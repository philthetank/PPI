//
//  main.c
//  PPI project  2D projection images from 2 opposed scintillator arrays
//
//  Created by Jurgen and Philip Seidel on 9/9/14.
//  Copyright (c) 2014 Jurgen and Philip Seidel. All rights reserved.
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
#define PROCESS_2D_PROJECTION_FLAG 'P'

//************************************** PROTOTYPES **************************
void    create_measured_planogram( float *planogram, const int planogram_size, char listfile[], char basename[]);
void    create_normplanogram( float *planogram, const int planogram_size);
void    parse_command_line( int argc, const char *argv[], char listfile[], char basename[] );
void    parse_from_stdin( char listfile[], char basename[] );
void    terminate_with_error( char message[] );
void    write_i2_array( unsigned short array[], int nelements, char filename[] );
void    write_i4_array( unsigned int array[], int nelements, char filename[] );
void    write_r4_array( float array[], int nelements, char filename[] );
void    write_array( const void *[], int, int, char[] );
void    process_2Dprojection_image( char[], char[] );
int     get_crystalindex( int );
int     calculate_projection_midplane_index( int, int, double );
void    sino_coord( double, double, double, double, double, double, double*, double*, double*, double*);


//************************************** GLOBALS **************************
double  detector_distance = 40;      // distance between surfaces of LYSO arrays
double  cone_angle = 5.0;       // cone angle in degrees;
double  crystal_pitch = 1.6;     // in mm, distance between crystal centers
int		keV_min = 250;
int		keV_max = 700;
int     conventional_2D_projection = 0;
double  acq_time_minutes = 1200.0; // acquisition time in minutes
double  bin_width = 0.8; // == half of crystal pitch
double  angular_bin_width = 2.1;
const int radial_bins = 2 * DETECTOR_COLS - 1;
const int axial_bins = 2 * DETECTOR_ROWS - 1;
const int projection_size = radial_bins * axial_bins;
const int azimuth_count = 17;
const int segment_count = 39;
const int span = 3;

//************************************** MAIN **************************
int main( int argc, const char * argv[] )
{
    char listfile[255], basename[255];
    
    // get input data file name
    if ( argc > 1) {
        // use command line arguments
        parse_command_line( argc, argv, listfile, basename);
        printf("Command line parsed.\n");
    } else {
        // parse arguments from stdin // currently only gets file path
        parse_from_stdin( listfile, basename);
        printf("File path parsed.\n");
    }
    
    if(conventional_2D_projection) {
        process_2Dprojection_image(listfile, basename);
        printf("Finished with conventional 2D projection image.\n");
        return 0;
    }
    
    const int planogram_size = projection_size * azimuth_count * segment_count;
    
    // create a "normalization" planogram
    float *normplanogram;
    normplanogram = calloc(planogram_size, sizeof(float));
    create_normplanogram( normplanogram, planogram_size);
    
    // create the planogram from measured data in list mode file
    float *measured_planogram;
    measured_planogram = calloc(planogram_size, sizeof(float));
    create_measured_planogram(measured_planogram, planogram_size, listfile, basename);
    
    // initialize image space
    int xbins = radial_bins;
    int ybins = round(detector_distance/0.8);
    int zbins = axial_bins;
    int image_size = xbins * ybins * zbins;
    float *image;
    image = calloc( image_size, sizeof(float));
    
    printf("Wrote images to disk.\n");
    printf("Finished with PLATO code.\n");
    return 0;
}


//************************************** create_normplanogram **************************
void create_measured_planogram( float *planogram, const int planogram_size, char listfile[], char basename[])
{   // we reading list mode data file on an event-by-event basis and placing events into projections
    // NOTE: row and columns have been switched compared to previous 2D-Projection mode version
    
    FILE *fp = fopen( listfile, "rb" );
    if ( !fp ) {
        char message[255];
        sprintf( message, "Error opening listfile %s\n", listfile );
        terminate_with_error(message);
    }
    
    const int acqtimemilliseconds = acq_time_minutes * 60000;
    int time, crystalindex1, crystalindex2, eventblock[EVENTBLOCK_SIZE];
    
    unsigned int rawimageD1[DETECTOR_ROWS * DETECTOR_COLS] = {0};
    unsigned int rawimageD2[DETECTOR_ROWS * DETECTOR_COLS] = {0};
    
    double y1 = (detector_distance * (-0.5));
    double y2 = (detector_distance * 0.5);
    double xcenter = (DETECTOR_COLS - 1) * 0.5;
    double zcenter = (DETECTOR_ROWS - 1) * 0.5;
    double r, phi, z, theta;
    double xcenterbins = xcenter * crystal_pitch /bin_width;
    double zcenterbins = zcenter * crystal_pitch /bin_width;
    
    int  event_count=0;
    int  energy1, energy2;
    
    while ( fread(&eventblock, sizeof(int), EVENT_ELEMENT_COUNT, fp) == EVENT_ELEMENT_COUNT ) {
        // get time, check if valid, skip special case
        time = eventblock[0];
        if ( time > acqtimemilliseconds ) {
            break;
        } else if ( time == -1 ) {
            continue;
        }
        
        // energy test first. both energies must be in window [keVmin, keVmax]
        energy1 = eventblock[1]>>16;
        if( energy1 < keV_min || energy1 > keV_max) continue;
        energy2 = eventblock[2]>>16;
        if( energy2 < keV_min || energy2 > keV_max) continue;
        
        crystalindex1 = get_crystalindex( eventblock[1] );
        crystalindex2 = get_crystalindex( eventblock[2] );
        rawimageD1[crystalindex1]++;
        rawimageD2[crystalindex2]++;
        
        int col1 = crystalindex1 / DETECTOR_ROWS;
        int col2 = crystalindex2 / DETECTOR_ROWS;
        int row1 = crystalindex1 % DETECTOR_ROWS;
        int row2 = crystalindex2 % DETECTOR_ROWS;
        row2 = (DETECTOR_ROWS - 1) - row2;
        
        double x1 = (col1 - xcenter) * crystal_pitch;
        double x2 = (col2 - xcenter) * crystal_pitch;
        double z1 = (row1 - zcenter) * crystal_pitch;
        double z2 = (row2 - zcenter) * crystal_pitch;
        
        int azimuth_segment = round( (double) (col2 - col1) / (double) span );
        if (azimuth_segment < 0) { azimuth_segment += azimuth_count; }
        if ( azimuth_segment < 0 || azimuth_segment >= azimuth_count ) {
            printf("azimuth segment = %d\n", azimuth_segment);
            terminate_with_error("azimuth segment is out of bounds");
        }
        int tilt_segment = round( (double) (row2 - row1) / (double) span );
        if ( tilt_segment < 0 ) { tilt_segment += segment_count; }
        if ( tilt_segment < 0 || tilt_segment >= segment_count ) {
            printf("Tilt segment = %d\n", tilt_segment);
            terminate_with_error("tilt segment is out of bounds");
        }
        
        sino_coord(x1, x2, y1, y2, z1, z2, &r, &phi, &theta, &z);
        r += xcenterbins;
        z += zcenterbins;
        // implementing bilinear interpolation here
        int   i = floor(r);
        int   j = floor(z);
        float fx = r-i;
        float fy = z-j;
        
        int index = j * radial_bins + i;
        index += azimuth_segment * projection_size;
        index += tilt_segment * azimuth_count * projection_size;
        if (index < planogram_size && index >= 0){
            planogram[index]   += (1.0-fx)*(1.0-fy);
            planogram[index+1] +=  fx*(1.0-fy);
            planogram[index+radial_bins]   +=  (1.0-fx)*fy;
            planogram[index+radial_bins+1] +=  fx*fy;
        } else {
            printf("Something's wrong with r %6.3lf and z %6.3lf\n", r, round(z));
        }
        event_count++;
        if( (event_count % 1000000)==0 ) printf("%10d events processed. scan time = %7.1lf minutes\n", event_count, (float) time/60000.0);
        
    }
    
    printf("%10d events processed. scan time = %7.1lf minutes\n", event_count, (float) time/60000.0);
    fclose( fp );
    
    char  filename[255];
    sprintf(filename,"%s_Det1_raw.img", basename);
    write_i4_array( &rawimageD1[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    sprintf(filename,"%s_Det2_raw.img", basename);
    write_i4_array( &rawimageD2[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    
    sprintf(filename,"%s_measured_projections.img", basename);
    write_r4_array( &planogram[0], planogram_size, filename);
    
}

//************************************** create_normplanogram **************************
void create_normplanogram( float *normplanogram, const int planogram_size)
{   //we are assigning q count to each line-of-response and bin it into the planogram
    // this can be used in lieu of a true normalization correction
    
    char    filename[255];
    double y1 = (detector_distance * (-0.5));
    double y2 = (detector_distance * 0.5);
    double xcenter = (DETECTOR_COLS - 1) * 0.5;
    double zcenter = (DETECTOR_ROWS - 1) * 0.5;
    double r, phi, z, theta;
    double xcenterbins = xcenter * crystal_pitch /bin_width;
    double zcenterbins = zcenter * crystal_pitch /bin_width;
    
    for (int row1 = 0; row1 < DETECTOR_ROWS; row1++) {
        double z1 = (row1 - zcenter) * crystal_pitch;
        for (int row2 = 0; row2 < DETECTOR_ROWS; row2++) {
            double z2 = (row2 - zcenter) * crystal_pitch;
            int tilt_segment = round( (double) (row2 - row1) / (double) span );
            if ( tilt_segment < 0 ) { tilt_segment += segment_count; }
            if ( tilt_segment < 0 || tilt_segment >= segment_count ) {
                printf("Tilt segment = %d\n", tilt_segment);
                terminate_with_error("tilt segment is out of bounds");
            }
            
            for (int col1 = 0; col1 < DETECTOR_COLS; col1++) {
                double x1 = (col1 - xcenter) * crystal_pitch;
                for (int col2 = 0; col2 < DETECTOR_COLS; col2++) {
                    int azimuth_segment = round( (double) (col2 - col1) / (double) span );
                    if (azimuth_segment < 0) { azimuth_segment += azimuth_count; }
                    if ( azimuth_segment < 0 || azimuth_segment >= azimuth_count ) {
                        printf("azimuth segment = %d\n", azimuth_segment);
                        terminate_with_error("azimuth segment is out of bounds");
                    }
                    double x2 = (col2 - xcenter) * crystal_pitch;
                    
                    sino_coord(x1, x2, y1, y2, z1, z2, &r, &phi, &theta, &z);
                    r += xcenterbins;
                    z += zcenterbins;
                    // implementing bilinear interpolation here
                    int  i = floor(r);
                    int  j = floor(z);
                    float fx = r-i;
                    float fy = z-j;
                    int index = j * radial_bins + i;
                    index += (tilt_segment * azimuth_count + azimuth_segment) * projection_size;
                    //                    index += azimuth_segment * projection_size;
                    //                    index += tilt_segment * azimuth_count * projection_size;
                    if (index < planogram_size && index >= 0){
                        normplanogram[index]   += (1.0-fx)*(1.0-fy);
                        normplanogram[index+1] +=  fx*(1.0-fy);
                        normplanogram[index+radial_bins]   +=  (1.0-fx)*fy;
                        normplanogram[index+radial_bins+1] +=  fx*fy;
                    } else {
                        printf("Somethings wrong with r %6.3lf and z %6.3lf\n", r, round(z));
                    }
                }
            }
        }
    }
    printf("norm planogram created.\n");
    sprintf(filename, "norm_planogram.img");
    write_r4_array(&normplanogram[0], planogram_size, filename);
}


//************************************** process_projection_image **************************
void process_2Dprojection_image(char listfile[], char basename[])
{
    // open file
    FILE *fp = fopen( listfile, "rb" );
    if ( !fp ) {
        char message[255];
        sprintf( message, "Error opening listfile %s\n", listfile );
        terminate_with_error(message);
    }
    
    unsigned int rawimageD1[DETECTOR_ROWS * DETECTOR_COLS] = {0};
    unsigned int rawimageD2[DETECTOR_ROWS * DETECTOR_COLS] = {0};
    unsigned int proj_image[PROJECTOR_ROWS * PROJECTOR_COLS] = {0};
    
    const int acqtimemilliseconds = acq_time_minutes * 60000;
    int     time, crystalindex1, crystalindex2;
    int     projectionindex, eventblock[EVENTBLOCK_SIZE];
    double acceptance_radius = (detector_distance/crystal_pitch)*tan(M_PI*cone_angle/180.0);
    double acceptance_radius2 = acceptance_radius*acceptance_radius;
    printf("cone_angle = %6.1lf degrees  distance = %8.1lf mm\n", cone_angle, detector_distance);
    printf("acceptance radius = %8.2lf  crystal units\n", acceptance_radius);
    
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
    
    char filename[255];
    sprintf(filename,"%s_Det1_raw.img", basename);
    write_i4_array( &rawimageD1[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    sprintf(filename,"%s_Det2_raw.img", basename);
    write_i4_array( &rawimageD2[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    sprintf(filename,"%s_projection.img", basename);
    write_i4_array( &proj_image[0], PROJECTOR_ROWS*PROJECTOR_COLS, filename);
}


//************************************** sino_coordinates **************************
void sino_coord( double x1, double x2, double y1, double y2, double z1, double z2, double* r, double* phi, double* theta, double *z)
{
    *phi = atan2(x2 - x1, y2 - y1);           // azimuth angle ( in rads )
    *r = (-y2 * sin(*phi) + x2 * cos(*phi));  // radial distance (in mm)
    //    *phi *= (180.0 / M_PI) / angular_bin_width;   // converts to angular bin number  (not needed currently)
    
    double temp = (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1);
    *theta = atan((z2 - z1) / sqrt(temp));     // tilt (or polar) angle
    
    *z =sqrt(y1 * y1 + x1 * x1 - (*r)*(*r)) * sin(*theta);  // axial location (in mm)
    *z += z1 * cos(*theta);
    *z /= bin_width;
    *r /= bin_width;
}


//************************************** calculate_projection_midplane_index **************************
int calculate_projection_midplane_index( int crystalindex1, int crystalindex2, double acceptance_radius2 )
{
    int col1 = crystalindex1 / DETECTOR_ROWS;
    int col2 = crystalindex2 / DETECTOR_ROWS;
    
    int row1 = crystalindex1 % DETECTOR_ROWS;
    int row2 = crystalindex2 % DETECTOR_ROWS;
    row2 = (DETECTOR_ROWS - 1) - row2;
    
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
                case PROCESS_2D_PROJECTION_FLAG:
                    printf("creating conventional 2D projection image.\n");
                    conventional_2D_projection = 1;
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

//************************************** write_r4_array ********************
void write_r4_array(float array[], int nelements, char filename[])
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
