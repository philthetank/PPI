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

#define INDEX_MASK  0x0000FFFF
#define ENERGY_MASK 0xFFFF0000

#define EVENT_ELEMENT_COUNT 3

#define ERROR_FILE_NAME "PPI_error_message.txt"

//************************************** PROTOTYPES **************************
void    parse_command_line( int argc, const char *argv[], char listfile[], char basename[] );
void    parse_from_stdin( char listfile[], char basename[] );
void    terminate_with_error( char message[] );
void    write_i2_array( unsigned short array[], int nelements, char filename[] );
void    write_i4_array( unsigned int array[], int nelements, char filename[] );
void    process_projection_image( char[], unsigned int[], unsigned int[], unsigned int[] );
int     get_crystalindex( int );
int     calculate_projection_midplane_index( int, int );

//************************************** GLOBALS **************************
double  D = 228.0;      // distance between surfaces of LYSO arrays
double  CA = 3.5;       // cone angle in degrees;
int		keVmin = 400;
int		keVmax = 650;
double  acqtime= 1200.0; // acquisition time in minutes

//************************************** MAIN **************************
int main( int argc, const char * argv[] )
{
    char listfile[255], filename[255], basename[255];

    unsigned int rawimageD1[DETECTOR_ROWS * DETECTOR_COLS] = {0};
    unsigned int rawimageD2[DETECTOR_ROWS * DETECTOR_COLS] = {0};
    unsigned int proj_image[PROJECTOR_ROWS * PROJECTOR_COLS] = {0};
    
    // get input data file name
    if ( argc > 1) {
        // use command line arguments
        parse_command_line( argc, argv, listfile, basename);
        printf("Command line parsed.\n");
    } else {
        // parse arguments from stdin // currently only gets file path
        parse_from_stdin( listfile, basename );
        printf("File path parsed.\n");
    }
    
    process_projection_image( listfile, rawimageD1, rawimageD2, proj_image );
    printf("Processed images.\n");
    
    sprintf(filename,"%s_Det1_raw.img", basename);
    write_i4_array( &rawimageD1[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    sprintf(filename,"%s_Det2_raw.img", basename);
    write_i4_array( &rawimageD2[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    sprintf(filename,"%s_projection.img", basename);
    write_i4_array( &proj_image[0], PROJECTOR_ROWS*PROJECTOR_COLS, filename);
    printf("Wrote images to disk.\n");
    
    printf("Finished with PPI.\n");
    return 0;
}

//************************************** process_projection_image **************************
void process_projection_image(char listfile[],
                              unsigned int rawimageD1[], unsigned int rawimageD2[],
                              unsigned int proj_image[] )
{
    // open file
    FILE *fp = fopen( listfile, "rb" );
    if ( !fp ) {
        char message[255];
        sprintf( message, "Error opening listfile %s\n", listfile );
        terminate_with_error(message);
    }
    
    const int acqtimemilliseconds = acqtime * 60000;
    int time, crystalindex1, crystalindex2, projectionindex, eventblock[EVENTBLOCK_SIZE];
    
    while ( fread(&eventblock, sizeof(int), EVENT_ELEMENT_COUNT, fp) == EVENT_ELEMENT_COUNT ) {
        // get time, check if valid, skip special case
        time = eventblock[0];
        if ( time > acqtimemilliseconds ) {
            break;
        } else if ( time == -1 ) {
            continue;
        }
        
        crystalindex1 = get_crystalindex( eventblock[1] );
        crystalindex2 = get_crystalindex( eventblock[2] );
        
        rawimageD1[crystalindex1]++;
        rawimageD2[crystalindex2]++;
        
        projectionindex = calculate_projection_midplane_index( crystalindex1, crystalindex2 );
        proj_image[projectionindex]++;
    }
    
    fclose( fp );
}

//************************************** calculate_projection_midplane_index **************************
int calculate_projection_midplane_index( int crystalindex1, int crystalindex2 )
{
    int row1 = crystalindex1 / DETECTOR_COLS;
    int row2 = crystalindex2 / DETECTOR_COLS;
    
    int col1 = crystalindex1 % DETECTOR_COLS;
    int col2 = crystalindex2 % DETECTOR_COLS;
    col2 = (DETECTOR_COLS - 1) - col2;
    
    int midplanerow = row1 + row2;
    int midplanecol = col1 + col2;
    
    int projectionindex = (midplanerow * PROJECTOR_COLS) + midplanecol;
    
    if ( projectionindex > MAX_INDEX ) {
        char message[255];
        sprintf(message, "Projection Index is out of bounds. Should be less than %d, but is %d\n", MAX_INDEX, projectionindex);
        terminate_with_error(message);
    }
    
    return projectionindex;
}

//************************************** get_crystalindex **************************
int get_crystalindex( int eventblock )
{
    int crystalindex = eventblock & INDEX_MASK;
    if ( crystalindex > MAX_INDEX ) {
        char message[255];
        sprintf(message, "Crystal Index is out of bounds. Should be less than %d, but is %d\n", MAX_INDEX, crystalindex);
        terminate_with_error(message);
    }
    return crystalindex;
}

//************************************** parse_command_line **************************
void parse_command_line(int argc, const char *argv[], char listfile[], char basename[])
{
    // command line must specify list file name first and it must be followed
    //    by the other parameters
    // e.g.  ./listproc LYSOtest.bin -a400 -b650
    
    if(argc <= 1) {
        terminate_with_error("Not enough arguments to find filename!");
    }
    
    if( sscanf( &argv[1][0], "%[^\t\n]", &listfile[0]) <= 0 ) {
        terminate_with_error("No luck parsing command line for filename!");
    }
    
    printf("Filename = %s\n", listfile);
    
    // basename = list file name without extension
    strcpy(basename, listfile);
    strtok(basename,".");
    
    puts(basename);
    puts(listfile);
    
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
                {
                    char message[255];
                    sprintf(message, "parse_command_line found unknown argument(s)!  arg =%c\n", argv[i][1]);
                    terminate_with_error(message);
                }
            }
        }
    }
}

//************************************** parse_from_stdin **************************
void parse_from_stdin( char listfile[], char basename[] )
{
    // get arguments interactively
    char buf[255];
    printf( "Enter listfile path:\n" );
    fgets( buf, 255, stdin );
    
    if ( sscanf( &buf[0], "%[^\t\n]", &listfile[0]) > 0 ) {
        strcpy( basename, listfile );
        strtok( basename,"." );
        
        puts( basename );
        puts( listfile );
    } else {
        terminate_with_error( "unable to get listfile path at runtime.\n" );
    }
}

//************************************** terminate with error **************************
void terminate_with_error(char message[])
{
    puts("Terminating with error: ");
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
    if ( !fp ) {
        char message[255];
        sprintf(message,"Error opening file %s\n", filename);
        terminate_with_error(message);
    }
    
    size_t n = fwrite( &array[0], 2, nelements, fp);
    if( n != nelements )  {
        char message[255];
        sprintf(message,"Write error! Only %zd words written instead of %d, file = %s\n", n, nelements, filename);
        terminate_with_error(message);
    }
    
    printf ("File %s was written to disk.\n", filename);
    fclose(fp);
}

//************************************** write_i4_array ********************
void write_i4_array(unsigned int array[], int nelements, char filename[])
{
    FILE *fp = fopen( filename,"wb");
    if ( !fp ) {
        char message[255];
        sprintf(message,"Error opening file %s\n", filename);
        terminate_with_error(message);
    }
    
    size_t n = fwrite( &array[0], 4, nelements, fp);
    if( n != nelements ) {
        char message[255];
        sprintf(message,"Write error! Only %zd words written instead of %d, file = %s\n", n, nelements, filename);
        terminate_with_error(message);
    }
    
    printf ("File %s was written to disk.\n", filename);
    fclose(fp);
}