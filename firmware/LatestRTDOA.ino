#include <Audio.h>

AudioInputI2SQuad i2s1;
AudioRecordQueue  queue1;   // Mic 1 (pin 8, L)
AudioRecordQueue  queue2;   //Mic 2 (pin 8, R)
AudioRecordQueue  queue3;   // Mic 3 (pin 6, L)
AudioConnection patchCord1(i2s1, 0, queue1, 0);
AudioConnection patchCord2(i2s1, 1, queue2, 0);
AudioConnection patchCord3(i2s1, 2, queue3, 0);


const int SAMPLES = 128;
const float TRIGGER = 200.0f;
const char     MAGIC[]       = "START";
const uint16_t RECORD_BLOCKS = 1034;                        // 1034 x 128 samples = 3.0 s , Sample rate = 44117Hz, 1block = 128 samples
const uint16_t BLOCKS_A      = 600;                         // first blocks live in RAM2
const uint16_t BLOCKS_B      = RECORD_BLOCKS - BLOCKS_A;    // the rest live in RAM1

// one block = mic1[128], mic2[128], mic3[128] back to back
DMAMEM int16_t storeA[BLOCKS_A][3][128];    // 460,800 bytes
int16_t        storeB[BLOCKS_B][3][128];    // 333,312 bytes

bool run_DOA = false;

inline int16_t (*blockPtr(uint16_t n))[128] {
    return (n < BLOCKS_A) ? storeA[n] : storeB[n - BLOCKS_A];
}

void armTrigger() {
    queue1.clear(); queue2.clear(); queue3.clear();
    queue1.begin(); queue2.begin(); queue3.begin();
    run_DOA = false;
}

void recordToRam() {
    digitalWrite(12, HIGH);                        // LED on = recording
    uint16_t n = 1;
    while (n < RECORD_BLOCKS) {
        if (queue1.available() && queue2.available() && queue3.available()) {
            int16_t (*dst)[128] = blockPtr(n);
            memcpy(dst[0], queue1.readBuffer(), 256); queue1.freeBuffer();
            memcpy(dst[1], queue2.readBuffer(), 256); queue2.freeBuffer();
            memcpy(dst[2], queue3.readBuffer(), 256); queue3.freeBuffer();
            n++;
        }
    }
    queue1.end(); queue2.end(); queue3.end();
    queue1.clear(); queue2.clear(); queue3.clear();
    digitalWrite(12, LOW);
}


void sendAll(const uint8_t* p, size_t n) {
    while (n > 0) {
        size_t w = Serial.write(p, n);
        p += w;
        n -= w;
        if (w == 0) yield();               // buffer full, let USB drain
    }
}

void sendFromRam() {
    uint8_t pkt[5 + 4 + 4 + 768];                           // START + counter + total + samples = 781 bytes
    const uint32_t total = RECORD_BLOCKS;
    memcpy(pkt, MAGIC, 5);
    memcpy(pkt + 9, &total, 4);
    for (uint32_t n = 0; n < RECORD_BLOCKS; n++) {
        memcpy(pkt + 5, &n, 4);
        memcpy(pkt + 13, blockPtr(n), 768);
        sendAll(pkt, sizeof(pkt));
    }
    Serial.send_now();
}

void setup() {
  Serial.begin(115200);
  pinMode(12, OUTPUT);

  AudioMemory(60);

  armTrigger();
}

void loop() {
  if (!run_DOA) {
        if (queue1.available() && queue2.available() && queue3.available()) {
            int16_t *b1 = queue1.readBuffer();
            int16_t *b2 = queue2.readBuffer();
            int16_t *b3 = queue3.readBuffer();

            float m1 = 0, m2 = 0, m3 = 0;
            for (int i = 0; i < SAMPLES; i++) { m1 += b1[i]; m2 += b2[i]; m3 += b3[i]; }
            m1 /= SAMPLES; m2 /= SAMPLES; m3 /= SAMPLES;

            bool trig = false;
            for (int i = 0; i < SAMPLES; i++) {
                if (fabsf(b1[i] - m1) > TRIGGER || fabsf(b2[i] - m2) > TRIGGER || fabsf(b3[i] - m3) > TRIGGER) { 
                  trig = true; 
                  break; 
                }
            }

            if (trig) {                                     // keep the triggering block as block 0
                int16_t (*dst)[128] = blockPtr(0);
                memcpy(dst[0], b1, 256);
                memcpy(dst[1], b2, 256);
                memcpy(dst[2], b3, 256);
                run_DOA = true;
            }
            queue1.freeBuffer(); queue2.freeBuffer(); queue3.freeBuffer();   // free BEFORE any clear()/end()
        }
    } else {
        recordToRam();      // ~3 s, nothing is sent while recording
        sendFromRam();      // then stream the finished recording
        delay(5000);         // let the sound die away before re-arming
        armTrigger();
    }
}
