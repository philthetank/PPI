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
#define DETECTOR_ROWS (59)
#define DETECTOR_COLS (26)

#define PROJECTOR_ROWS ((2 * DETECTOR_ROWS) - 1)
#define PROJECTOR_COLS ((2 * DETECTOR_COLS) - 1)

#define EVENTBLOCK_SIZE (3 * sizeof(int))

#define DETECTOR_MAX_INDEX ((DETECTOR_ROWS * DETECTOR_COLS) - 1)
#define PROJECTOR_MAX_INDEX ((PROJECTOR_ROWS * PROJECTOR_COLS) - 1)

#define INDEX_MASK  0x0000FFFF

#define EVENT_ELEMENT_COUNT 3

#define ERROR_FILE_NAME "PPI_error_message.txt"

#define KEV_MIN_FLAG 'a'
#define KEV_MAX_FLAG 'b'
#define CONE_ANGLE_FLAG 'C'
#define DETECTOR_DISTANCE_FLAG 'D'
#define ACQ_TIME_MINUTES_FLAG 'm'

//************************************** PROTOTYPES **************************
void    parse_command_line( int argc, const char *argv[], char listfile[], char basename[] );
void    parse_from_stdin( char listfile[], char basename[] );
void    terminate_with_error( char message[] );
void    write_i2_array( unsigned short array[], int nelements, char filename[] );
void    write_i4_array( unsigned int array[], int nelements, char filename[] );
void    write_array( const void *[], int, int, char[] );
void    process_projection_image( char[], unsigned int[], unsigned int[], unsigned int[], double );
int     get_crystalindex( int );
int     calculate_projection_midplane_index( int, int, double );
void    sino_coord( double, double, double, double, double*, double*);

//************************************** GLOBALS **************************
double  detector_distance = 40;      // distance between surfaces of LYSO arrays
double  cone_angle = 3.5;       // cone angle in degrees;
double  crystal_pitch = 1.6;     // in mm, distance between crystal centers
int		keV_min = 400;
int		keV_max = 650;
double  acq_time_minutes = 1200.0; // acquisition time in minutes
double  bin_width = 0.8; // == half of crystal pitch
double  angular_bin_width = 2.1;

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
    
    /* not currently needed for sinogram. */
    
    //    double acceptance_radius = (detector_distance/crystal_pitch)*tan(M_PI*cone_angle/180.0);
    //    double acc_radius2 = acceptance_radius*acceptance_radius;
    //    printf("cone_angle = %6.1lf degrees  distance = %8.1lf mm\n", cone_angle, detector_distance);
    //    printf("acceptance radius = %8.2lf  crystal units\n", acceptance_radius);
    
    //    process_projection_image( listfile, rawimageD1, rawimageD2, proj_image, acc_radius2 );
    //    printf("Processed images.\n");
    
    double y1 = (detector_distance * (-0.5));
    double y2 = (detector_distance * 0.5);
    double center = (DETECTOR_COLS - 1) * 0.5;
    
    const int radial_bins = 51;
    const int angular_bins = 50;
    const int slices = 2 * DETECTOR_ROWS - 1;
    const int sino_size = radial_bins * angular_bins * slices;
    const int segment_count = 37;
    
    unsigned int sinogram[sino_size][segment_count] = {0};
    
    for (int row1 = 0; row1 < DETECTOR_ROWS; row1++) {
        for (int row2 = 0; row2 < DETECTOR_ROWS; row2++) {
            if ( abs(row1 - row2) > 1) { continue; }
            int slice_number = row1 + row2;
            for (int col1 = 0; col1 < DETECTOR_COLS; col1++) {
                double x1 = (col1 - center) * crystal_pitch;
                for (int col2 = 0; col2 < DETECTOR_COLS; col2++) {
                    double r, phi;
                    double x2 = (col2 - center) * crystal_pitch;
                    
                    sino_coord(x1, x2, y1, y2, &r, &phi);
                    phi += 25;
                    int radial_coord = round(r + angular_bins / 2);
                    int index = round(phi) * radial_bins + radial_coord;
                    index += radial_bins * angular_bins * slice_number;
                    
                    if (index < sino_size && index >= 0){
                        sinogram[index][0]++;
                        printf("These seem okay: r %6.3lf and phi %6.3lf\n", r, phi);
                    } else {
                        printf("Somethings wrong with r %6.3lf and phi %6.3lf\n", r, phi);
                    }
                    
                }
                printf("Column %d done\n", col1);
            }
        }
    }
    
    sprintf(filename, "%s_sinogram.img", basename);
    write_i4_array(&sinogram[0], sino_size, filename);
    
    //    sprintf(filename,"%s_Det1_raw.img", basename);
    //    write_i4_array( &rawimageD1[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    //    sprintf(filename,"%s_Det2_raw.img", basename);
    //    write_i4_array( &rawimageD2[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    //    sprintf(filename,"%s_projection.img", basename);
    //    write_i4_array( &proj_image[0], PROJECTOR_ROWS*PROJECTOR_COLS, filename);
    printf("Wrote images to disk.\n");
    
    printf("Finished with PPI.\n");
    return 0;
}

//************************************** process_projection_image **************************
void process_projection_image(char listfile[],
                              unsigned int rawimageD1[], unsigned int rawimageD2[],
                              unsigned int proj_image[],
                              double acceptance_radius2)
{
    // open file
    FILE *fp = fopen( listfile, "rb" );
    if ( !fp ) {
        char message[255];
        sprintf( message, "Error opening listfile %s\n", listfile );
        terminate_with_error(message);
    }
    
    const int acqtimemilliseconds = acq_time_minutes * 60000;
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
        
        projectionindex = calculate_projection_midplane_index( crystalindex1, crystalindex2, acceptance_radius2 );
        if ( projectionindex >= 0) {
            proj_image[projectionindex]++;
        }
    }
    
    fclose( fp );
}


void sino_coord( double x1, double x2, double y1, double y2, double* r, double* phi)
{
    *phi = atan2(x2 - x1, y2 - y1); // gets angle ( in rads )
    *r = (-y2 * sin(*phi) + x2 * cos(*phi)) / bin_width; // get distance in bin width
    *phi *= (180.0 / M_PI) / angular_bin_width; // converts to angular bin number
}

//************************************** calculate_projection_midplane_index **************************
int calculate_projection_midplane_index( int crystalindex1, int crystalindex2, double acceptance_radius2 )
{
    int row1 = crystalindex1 / DETECTOR_COLS;
    int row2 = crystalindex2 / DETECTOR_COLS;
    
    int col1 = crystalindex1 % DETECTOR_COLS;
    int col2 = crystalindex2 % DETECTOR_COLS;
    col2 = (DETECTOR_COLS - 1) - col2;
    
    int midplanerow = row1 + row2;
    int midplanecol = col1 + col2;
    
    int distance2 = ((row1 - row2) * (row1 - row2)) + ((col1 - col2) * (col1 - col2));
    int projectionindex = (distance2 > acceptance_radius2) ? -1 : (midplanerow * PROJECTOR_COLS) + midplanecol;
    
    if ( projectionindex > PROJECTOR_MAX_INDEX ) {
        char message[255];
        sprintf(message, "Projection Index is out of bounds. Should be less than %d, but is %d\n", PROJECTOR_MAX_INDEX, projectionindex);
        terminate_with_error(message);
    }
    
    return projectionindex;
}

//************************************** get_crystalindex **************************
int get_crystalindex( int eventblock )
{
    int crystalindex = eventblock & INDEX_MASK;
    if ( crystalindex > DETECTOR_MAX_INDEX ) {
        char message[255];
        sprintf(message, "Crystal Index is out of bounds. Should be less than %d, but is %d\n", DETECTOR_MAX_INDEX, crystalindex);
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
    char * extension = strrchr(listfile,'.');
    if(extension == NULL)
        strcpy( basename, listfile);
    else {
        // len is the length of the basename
        size_t len = strlen(listfile) - strlen(extension);
        strncpy( basename, listfile, len);
    }
    printf("basename = %s  extension = %s\n", basename, extension);
    
    for (int i=2; i<argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case KEV_MIN_FLAG:								// lower threshold (must be above a, in keV)
                    keV_min = atoi(&argv[i][2]);
                    printf("lower threshold (keV) = %3d\n", keV_min);
                    break;
                case KEV_MAX_FLAG:                               // upper threshold (must be below b, in keV)
                    keV_max = atoi(&argv[i][2]);
                    printf("upper threshold (keV) = %3d\n", keV_max);
                    break;
                case CONE_ANGLE_FLAG:                               // cone angle (in degrees)
                    sscanf(&argv[i][2],"%lf", &cone_angle);
                    printf("cone angle =%6.2lf degrees\n", cone_angle);
                    break;
                case DETECTOR_DISTANCE_FLAG:                               // distance (in millimeters)
                    sscanf(&argv[i][2],"%lf", &detector_distance);
                    printf("distance between crystal surfaces =%8.2lf mm\n", detector_distance);
                    break;
                case ACQ_TIME_MINUTES_FLAG:                               // acquisition time in minutes
                    sscanf(&argv[i][2],"%lf", &acq_time_minutes);
                    printf("acquisition time = %6.1lf minutes\n", acq_time_minutes);
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
        // basename = list file name without extension
        char * extension = strrchr(listfile,'.');
        if(extension == NULL)
            strcpy( basename, listfile);
        else {
            // len is the length of the basename
            size_t len = strlen(listfile) - strlen(extension);
            strncpy( basename, listfile, len);
        }
        printf("basename = %s  extension = %s\n", basename, extension);
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

//************************************** write_array ********************
void write_array( const void *array[], int size, int count, char filename[] )
{
    FILE *fp = fopen( filename,"wb");
    if ( !fp ) {
        char message[255];
        sprintf(message,"Error opening file %s\n", filename);
        terminate_with_error(message);
    }
    
    size_t n = fwrite( &array[0], size, count, fp);
    if( n != count ) {
        char message[255];
        sprintf(message,"Write error! Only %zd words written instead of %d, file = %s\n", n, count, filename);
        terminate_with_error(message);
    }
    
    printf ("File %s was written to disk.\n", filename);
    fclose(fp);
}
