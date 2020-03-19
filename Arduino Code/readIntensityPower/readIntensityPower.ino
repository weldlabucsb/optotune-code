//Author: Max
//Date: 7/25/2019

//This is some simple code to retrive analog values from the arduino
//monitoring the voltages on the photodetectors


int powerPin = A1;
int intensityPin = A0;
int power = 0;
int intensity = 0;


void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  Serial.flush();
}

void serialEvent() { 
  char inChar=(char)Serial.read();
  if (inChar=='\n' || inChar=='\r'){
    readValue();
  }
}


void readValue(){
  power = analogRead(powerPin);
  intensity = analogRead(intensityPin);
  Serial.print(power,DEC);
  Serial.write('\r');
  Serial.print(intensity,DEC);
  Serial.write('\r');
}
void loop() {
  // put your main code here, to run repeatedly:

}
