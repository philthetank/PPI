//
//  main.c
//  PPI project  2D projection images from 2 opposed scintillator arrays
//
//  Created by Jurgen Seidel on 9/9/14.
//  Copyright (c) 2014 Jurgen Seidel. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//************************************** CONSTANTS **************************
#define DETECTOR_ROWS 26
#define DETECTOR_COLS 59

#define ERROR_FILE "PPI_error_message.txt"

//************************************** PROTOTYPES **************************
void    parse_command_line(int argc, const char *argv[], char listfile[], char basename[]);
void    terminate_with_error( char message[]);
void    write_i2_array( unsigned short array[], int nelements, char filename[]);


//************************************** GLOBALS **************************
double  D = 228.0;      // distance between surfaces of LYSO arrays
double  CA = 3.5;       // cone angle in degrees;
int		keVmin = 400;
int		keVmax = 650;
double  acqtime= 1200.0; // acquisition time in minutes


//************************************** MAIN **************************
int main(int argc, const char * argv[])
{
    char 	listfile[255], filename[255], basename[255];
    unsigned short rawimageD1[DETECTOR_ROWS*DETECTOR_COLS];
    unsigned short rawimageD2[DETECTOR_ROWS*DETECTOR_COLS];
    
    // get input data file name
    parse_command_line( argc, argv, listfile, basename);
    printf("command line parsed.\n");
    
    sprintf(filename,"%s_Det1_raw.img", basename);
    write_i2_array( &rawimageD1[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    sprintf(filename,"%s_Det2_raw.img", basename);
    write_i2_array( &rawimageD2[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    
    return 0;
}

//************************************** parse_command_line **************************
void parse_command_line(int argc, const char *argv[], char listfile[], char basename[])
{
    // command line must specify list file name first and it must be followed
    //    by the other parameters
    // e.g.  ./listproc LYSOtest.bin -a400 -b650
    
    if(argc > 1) {
        int nrd = sscanf( &argv[1][0], "%s", &listfile[0]);
        if(nrd > 0) {
            printf("Filename = %s\n", listfile);
            
            // basename = list file name without extension
            strcpy(basename, listfile);
            strtok(basename,".");
            
            puts(basename);
            puts(listfile);
            
            char message[255];
            for (int i=2; i<argc; i++) {
                if (argv[i][0] == '-') {
                    switch(argv[i][1]) {
                        case 'a':								// lower threshold (must be above a, in keV)
                            keVmin = atoi(&argv[i][2]);
                            printf("lower threshold (keV) = %3d\n", keVmin);
                            break;
                        case 'b':                               // upper threshold (must be below b, in keV)
                            keVmax = atoi(&argv[i][2]);
                            printf("upper threshold (keV) = %3d\n", keVmax);
                            break;
                        case 'C':                               // cone angle (in degrees)
                            sscanf(&argv[i][2],"%lf", &CA);
                            printf("cone angle =%6.2lf degrees\n", CA);
                            break;
                        case 'D':                               // distance (in millimeters)
                            sscanf(&argv[i][2],"%lf", &D);
                            printf("distance between crystal surfaces =%8.2lf mm\n", D);
                            break;
                        case 'm':                               // acquisition time in minutes
                            sscanf(&argv[i][2],"%lf", &acqtime);
                            printf("acquisition time = %6.1lf minutes\n", acqtime);
                            break;
                        default:
                            sprintf(message, "parse_command_line found unknown argument(s)!  arg =%c\n", argv[i][1]);
                            terminate_with_error(message);
                    }
                }
            }
        } else {
            terminate_with_error("No luck parsing command line for filename!");
        }
    } else {
        terminate_with_error("Not enough arguments to find filename!");
    }
}

//************************************** terminate with error **************************
void terminate_with_error(char message[])
{
    puts(" ");
    puts(message);
    
    FILE *fp = fopen(ERROR_FILE,"w");
    fprintf(fp, "%s\n",message);
    fclose(fp);
    exit(-1);
}

//************************************** write_i2_array ********************
void write_i2_array(unsigned short array[], int nelements, char filename[])
{
    FILE *fp = fopen( filename,"wb");
    char message[255];
    
    if (fp) {
        size_t n = fwrite( &array[0], 2, nelements, fp);
        if(n != nelements)  {
            sprintf(message,"Write error! Only %zd words written instead of %d, file = %s\n", n, nelements, filename);
            terminate_with_error(message);
        }
        printf ("File  %s  was written to disk.\n", filename);
    } else {
        sprintf(message,"Error opening file %s\n", filename);
        terminate_with_error(message);
    }
    
    fclose(fp);
}
