///////////Use information//////////////
// This code is written to generate discrete voltages to drive ETL and synchronised with camera////////
// All the components must be connected as discribed in the article//////
// In ETL set lower and upper current limit and set Analog input mode of operation///// 
// Enble external trigger of the camera//////
//upload in arduino.... starts data taking/////

#include <SPI.h>
byte address = 0x11;
int CS= 10;
int i=0;
int ex=30;//set camera exposure time
int d=40; // setcamera readout time
int k=10; //  set number of planes in one cycle
int n=100; // set number of cycles for cyclic scanning 

void setup()
{
  pinMode (CS, OUTPUT);
  SPI.begin();
  // adjust high and low resistance of potentiometer
  // adjust Highest Resistance .
   digitalPotWrite(0x00);
   //delay(1000);

      // adjust  wiper in the  Mid point  .
   digitalPotWrite(0x80);
   //delay(1000);

   // adjust Lowest Resistance .
   digitalPotWrite(0xFF);
   //delay(1000);
    ////  camera trigger pin/////////
   pinMode(3, OUTPUT);  //camera is conected to Pin-3 of arduino 
   digitalWrite(3, LOW); // set camera trigger line LOW

   
}

void loop()
{

  for (int j=1; j<=n; j++){
    
for (int i = 1; i <= k; i++){
  rt(i);

}

}
exit(0); //  stop upon completion
    
}

int digitalPotWrite(int value)
{
  digitalWrite(CS, LOW);
  SPI.transfer(address);
  SPI.transfer(value);
  digitalWrite(CS, HIGH);
}



void rt(int i)
{
digitalPotWrite(2.55*i); //255 stand mimimum resiatance and 0 stand maximum resiatance of digital potentio meter   
//here trigger////////
digitalWrite(3, HIGH);   // start expose
delay(ex);                // exposure camera
digitalWrite(3, LOW);     // 
delay(d);  //readtime plus processing time

  }
