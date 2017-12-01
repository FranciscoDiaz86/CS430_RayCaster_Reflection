B## Student Info
Name: Francisco Diaz
Email: fd86@nau.edu
Class: CS 430

## Usage
This program takes 4 commands: The length, width, the CVS file, and the output PPM file. This version of the raytracer is now recursive so that it can calculate reflections. When I ray is shot out of the camera it will hit an object then
a reflection ray will be calculated and the new reflection ray will be shot out and look for other objects until it reaches the recursive limit or the t value is infinity. The reflection rays return the color of the objects it hits and it blends it in with the
closest object hit and you get your reflection.


###Example command line input
./raycast 100 100 test.csv output.ppm

## Known Issues
The issue that I have with this current version of the raytracer is that for some reason only the plane gets the reflection. When I try to make the sphere the reflective object you do not get a reflection on the sphere. 