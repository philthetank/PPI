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
#define DETECTOR_ROWS (26)
#define DETECTOR_COLS (59)

#define PROJECTOR_ROWS ((2 * DETECTOR_ROWS) - 1)
#define PROJECTOR_COLS ((2 * DETECTOR_COLS) - 1)

#define EVENTBLOCK_SIZE (3 * sizeof(int))

#define MAX_INDEX ((DETECTOR_ROWS * DETECTOR_COLS) - 1)

#define ERROR_FILE_NAME "PPI_error_message.txt"

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
    char listfile[255], filename[255], basename[255];
    unsigned short rawimageD1[DETECTOR_ROWS * DETECTOR_COLS];
    unsigned short rawimageD2[DETECTOR_ROWS * DETECTOR_COLS];
    unsigned short proj_image[PROJECTOR_ROWS * PROJECTOR_COLS];
    
    // get input data file name
    parse_command_line( argc, argv, listfile, basename);
    printf("Command line parsed.\n");
    
    // open file
    FILE *fp = fopen( listfile, "rb");
    if (fp) {
        int eventblock[EVENTBLOCK_SIZE];
        const int acqtimeMinutes = acqtime * 60;
        char message[255];
        while ( fread(&eventblock, EVENTBLOCK_SIZE, 1, fp) >= 3 ) {
            int time = eventblock[0];
            if ( time <= acqtimeMinutes ) {
                int crystalindex1 = eventblock[1] & 0x00FF;
                int crystalindex2 = eventblock[2] & 0x00FF;
                
                if ( crystalindex1 > MAX_INDEX ) {
                    sprintf(message, "Crystal Index #1 is out of bounds. Should be less than %d, but is %d\n", MAX_INDEX, crystalindex1);
                    terminate_with_error(message);
                }
                
                if ( crystalindex2 > MAX_INDEX ) {
                    sprintf(message, "Crystal Index #2 is out of bounds. Should be less than %d, but is %d\n", MAX_INDEX, crystalindex2);
                    terminate_with_error(message);
                }
                
                rawimageD1[crystalindex1]++;
                rawimageD2[crystalindex2]++;
                
                int row1 = crystalindex1 / DETECTOR_COLS;
                int row2 = crystalindex2 / DETECTOR_COLS;
                row2 = (DETECTOR_ROWS - 1) - row2;
                
                int col1 = crystalindex1 % DETECTOR_COLS;
                int col2 = crystalindex2 % DETECTOR_COLS;
                
                int midplanerow = row1 + row2;
                int midplanecol = col1 + col2;
                
                int projectionindex = midplanerow * PROJECTOR_COLS + midplanecol;
                
                proj_image[projectionindex]++;
                
                //                    int energy1 = eventblock[1] & 0xFF00;
                //                    int energy2 = eventblock[2] & 0xFF00;
                
            } else {
                break;
            }
        }
        fclose( fp );
    } else {
        char message[255];
        sprintf(message,"Error opening listfile %s\n", filename);
        terminate_with_error(message);
    }
    
    printf("Processed image.\n");
    
    sprintf(filename,"%s_Det1_raw.img", basename);
    write_i2_array( &rawimageD1[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    sprintf(filename,"%s_Det2_raw.img", basename);
    write_i2_array( &rawimageD2[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    sprintf(filename,"%s_projection.img", basename);
    write_i2_array( &proj_image[0], PROJECTOR_ROWS*PROJECTOR_COLS, filename);
    
    printf("Wrote image to disk.\n");

    printf("Finished with PPI.\n");
    return 0;
}

//************************************** parse_command_line **************************
void parse_command_line(int argc, const char *argv[], char listfile[], char basename[])
{
    // command line must specify list file name first and it must be followed
    //    by the other parameters
    // e.g.  ./listproc LYSOtest.bin -a400 -b650
    
    if(argc > 1) {
        int nrd = sscanf( &argv[1][0], "%[^\t\n]", &listfile[0]);
        
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
    
    FILE *fp = fopen(ERROR_FILE_NAME,"w");
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
