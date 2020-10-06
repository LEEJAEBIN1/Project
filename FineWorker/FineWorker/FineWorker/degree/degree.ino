#include <SoftwareSerial.h>
#include "Wire.h"
#include "MPU9250.h"
#include "I2Cdev.h"

SoftwareSerial BTSerial(3,2);//rx = 2번 tx = 3번

#define pi 3.141592
#define RADIANS_TO_DEGREES 180/3.14159
#define fs 131.0;
MPU9250 mpu;
int a0 = 0;
int Degree = 0;
int State = 0;
void Func_walking();
void Check_degree();


int16_t ax,ay,az;
int16_t gx,gy,gz;
int16_t mx,my,mz;
float gyro_z;
float accel_x;
//자이로센서 바이어스값
float base_gx=0, base_gy=0, base_gz=0;
unsigned long pre_msec=0;
void calibrate()
{  
  int loop =10;
  for (int i=0;i<loop;i++)
  {
    mpu.getMotion9(&ax,&ay,&az,&gx,&gy,&gz,&mx,&my,&mz);
    base_gx += gx;
    base_gy += gy;
    base_gz += gz;
    delay(200);
  }
  base_gx /=loop;
  base_gy /=loop;
  base_gz /=loop;
}
 
 
void setup() {
  Wire.begin();
  BTSerial.begin(9600);
  Serial.begin(9600);  
  mpu.initialize();
  calibrate(); 
}
void Check_degree()
{
    float dt = (millis()-pre_msec)/1000.0;
    pre_msec = millis();
    mpu.getMotion9(&ax,&ay,&az,&gx,&gy,&gz,&mx,&my,&mz);
    accel_x = atan(ay/sqrt(pow(ax,2) + pow(az,2)))*RADIANS_TO_DEGREES;
    gyro_z = (gz-base_gz)/fs;
}
void loop() {
  Check_degree();
  Serial.print((int)accel_x);
  Serial.print("\t");
  Serial.println((int)gyro_z);
  BTSerial.write("12");
}
