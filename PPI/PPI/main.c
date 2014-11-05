//
//  main.c
//  PLATO project  3D image made with 2 opposed scintillator arrays
//
//  Created by Jurgen and Philip Seidel on 9/9/14.
//  Copyright (c) 2014 Jurgen and Philip Seidel. All rights reserved.
//
//  new: input can be either a list mode file full of event data (extension ".bin")
//                        or a projection data file  (file name ends with "measured_projections.img")
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

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
#define NORMALIZATION_FLAG 'N'
#define ITERATIONS_FLAG 'i'
#define TEST_FLAG 'T'

//#define	MIN(a,b)  ( ((a)<(b)) ? (a) : (b))
//#define	MAX(a,b)  ( ((a)>(b)) ? (a) : (b))


//************************************** PROTOTYPES **************************
void  back_project( float *image, short *image_mask, float *factors, float *weights, int xbins, int ybins, int zbins,
                   float *measured_planogram, float *est_planogram, int planogram_size);
void    calculate_normalization_correction( float *measured_planogram, int planogram_size);
void    create_measured_planogram( float *planogram, const int planogram_size, char listfile[], char basename[]);
void    create_normplanogram( float *planogram, const int planogram_size);
void    create_image_mask( short int *image_mask, int xbins_half, int ybins_half, int zbins_half, float *normplanogram, int planogram_size);
void    forward_project( float *image, short int *image_mask, int xbins, int ybins, int zbins, float *est_planogram, int planogram_size);
void    parse_command_line( int argc, const char *argv[], char listfile[], char basename[], char path[] );
void    parse_from_stdin( char listfile[], char basename[], char path[] );
void    terminate_with_error( char message[] );
void    read_r4_array(float array[], int nelements, char filename[]);
void    write_i2_array( short array[], int nelements, char filename[] );
void    write_i4_array( unsigned int array[], int nelements, char filename[] );
void    write_r4_array( float array[], int nelements, char filename[] );
void    write_array( const void *[], int, int, char[] );
void    process_2Dprojection_image( char[], char[] );
int     get_crystalindex( int );
int     calculate_projection_midplane_index( int, int, double );
void    sino_coord( double, double, double, double, double, double, double*, double*, double*, double*);
float   gaussian( double x, double FWHM);


//************************************** GLOBALS **************************
double  detector_distance = 40;      // distance between surfaces of LYSO arrays
double  cone_angle = 5.0;       // cone angle in degrees;
double  crystal_pitch = 1.6;     // in mm, distance between crystal centers
int		keV_min = 250;
int		keV_max = 700;
int     conventional_2D_projection = 0;
int     process_normalization = 0;
double  acq_time_minutes = 1200.0; // acquisition time in minutes
double  bin_width = 0.8; // == half of crystal pitch
double  angular_bin_width = 2.1;
const int radial_bins = 2 * DETECTOR_COLS - 1;
const int axial_bins = 2 * DETECTOR_ROWS - 1;
const int projection_size = 51*117;  //radial_bins * axial_bins;
const int azimuth_count = 17;  // could be as high as 51?
const int segment_count = 39;
const int span = 3;
double  voxel_size = 0.8;   // later change it to 0.4 mm
int     iterations = 50;
int     testing = 0;
int     is_listfile = 0;

//************************************** MAIN **************************
int main( int argc, const char * argv[] )
{
    char listfile[255], filename[255], basename[255], path[255];
    
    // get input data file name
    if ( argc > 1) {
        // use command line arguments
        parse_command_line( argc, argv, listfile, basename, path);
        printf("Command line parsed.\n");
    } else {
        // parse arguments from stdin // currently only gets file path
        parse_from_stdin( listfile, basename, path);
        printf("File path parsed.\n");
    }
    
    if(conventional_2D_projection) {
        process_2Dprojection_image(listfile, basename);
        printf("Finished with conventional 2D projection image.\n");
        return 0;
    }
    
    if(testing) {
        const int planogram_size = projection_size * 51 * segment_count;
        float *normplanogram = calloc(planogram_size, sizeof(float));
        create_normplanogram( normplanogram, planogram_size);
        return 0;
    }
    
    const int planogram_size = projection_size * azimuth_count * segment_count;
    float *normplanogram;
    float *measured_planogram = calloc(planogram_size, sizeof(float));

    if(is_listfile) {      // bin the measured data (from the list mode file) into planogram
        create_measured_planogram( measured_planogram, planogram_size, listfile, basename);
    }
    else {      // read input planogram from file
        read_r4_array( measured_planogram, planogram_size, listfile);
    }
    

    // create a "normalization" planogram or else read the normalization from file
    if(process_normalization == 1) {
        calculate_normalization_correction( measured_planogram, planogram_size);
        printf("finished with normalization procedure.\n");
        return 0;
    }
    else {   // read normalization file
        normplanogram = calloc(planogram_size, sizeof(float));
        char normalization_file[255];
        sprintf( normalization_file, "%snormalization_planogram.img", path);
        read_r4_array( normplanogram, planogram_size, normalization_file);
    }
    
    // set up image space
    int xbins = radial_bins;                    // radial_bins = 51  (40 mm)
    int ybins = round(0.9*detector_distance/0.8);   //  ybins  = 45  (36 mm) we are using 90% of the detector_distance
    int zbins = axial_bins;                     // axial_bins = 117  (93 mm)
    int image_size = xbins * ybins * zbins;
    float *image = calloc( image_size, sizeof(float));
    
    //initialize image
    double totalcounts = 0.0;
    for(int i = 0; i<planogram_size; i++) { totalcounts += measured_planogram[i]; }
    float counts_per_voxel = totalcounts/image_size;
    for(int i = 0; i<image_size; i++) { image[i] = counts_per_voxel; }
    
    // time variables
    struct 	timeval  tv_start;
    struct 	timezone tz_start;
    // start timer
    gettimeofday( &tv_start, &tz_start);
    double tstart = tv_start.tv_sec + 0.000001*tv_start.tv_usec;
    
    // create image mask that tells us which part of the image will show up in a given planogram
    int image_mask_size =  xbins*ybins*zbins;
    image_mask_size *= segment_count*azimuth_count;
    printf("Size of image mask = %d bytes\n", image_mask_size*2);
    short int *image_mask = calloc( image_mask_size, sizeof(short int));
    create_image_mask( image_mask, xbins, ybins, zbins, normplanogram, planogram_size);
    
    gettimeofday(&tv_start,&tz_start);
    double tstop = tv_start.tv_sec + 0.000001*tv_start.tv_usec;
    double elapsed_time = tstop - tstart;
    printf("create image mask took %10.2lf seconds.\n", elapsed_time);
    double tstart2 = tv_start.tv_sec + 0.000001*tv_start.tv_usec;
    
    float *est_planogram = calloc(planogram_size, sizeof(float));
    float *factors = calloc( image_size, sizeof(float));
    float *weights = calloc( image_size, sizeof(float));
    
    // start the cycle over iterations
    for(int i=0; i<iterations; i++) {
        forward_project( image, image_mask, xbins, ybins, zbins, est_planogram, planogram_size);

        sprintf(filename,"%sestimated projections iteration_%2.2d.raw", path, i+1);
        write_r4_array( est_planogram, planogram_size, filename);

        for( int j=0; j< planogram_size; j++) est_planogram[j] *= normplanogram[j];

        back_project( image, image_mask, factors, weights, xbins, ybins, zbins, measured_planogram, est_planogram, planogram_size);
        sprintf(filename,"%sestimated image iteration_%2.2d.raw", path, i+1);
        write_r4_array( image, xbins*ybins*zbins, filename);
    }
    
/*    sprintf(filename,"%sestimated projections iteration_%2.2d.raw", path, iterations);
    write_r4_array( est_planogram, planogram_size, filename);
    sprintf(filename,"%sestimated image iteration_%2.2d.raw", path, iterations);
    write_r4_array( image, xbins*ybins*zbins, filename);
*/
    gettimeofday(&tv_start,&tz_start);
    tstop = tv_start.tv_sec + 0.000001*tv_start.tv_usec;
    printf("\n%2d iteration(s) took %10.2lf seconds.\n", iterations, tstop-tstart2);
    printf("\ndone!\n");
    
    printf("Finished with PLATO code.\n\n");
    return 0;
}


//************************************** back projection *****************************************
void  back_project( float *image, short *image_mask, float *factors, float *weights, int xbins, int ybins, int zbins,
                   float *measured_planogram, float *est_planogram, int planogram_size)
{
    // takes the ratio of measured to estimated projections (point by point) and projects it back over the image space
    // ratio = measured_planogram[i]/est_planogram[i];  (if measured_planogram > 0; else 0)
    
    double xcenter = (DETECTOR_COLS - 1) * 0.5;
    double zcenter = (DETECTOR_ROWS - 1) * 0.5;
    double xcenterbins = xcenter * crystal_pitch /bin_width;
    double zcenterbins = zcenter * crystal_pitch /bin_width;
    double temp, theta, phi, x,y,z, r,s;
    double dx, dz, FWHM=1.5;        // FWHM in bins  = 1.5*0.8 mm = 1.2 mm
    float  wtx[radial_bins];
    
    // convert estimated planogram to ratio of measured to estimated.
    for( int i=0; i< planogram_size; i++) {
        if((est_planogram[i]) > 0.01) {
            est_planogram[i] = measured_planogram[i]/est_planogram[i];
        }
        else {
            est_planogram[i] = 0.0;
        }
    }
    
    for( int i=0; i< xbins*ybins*zbins; i++) factors[i]=weights[i]=0.0;
    
    // for all azimuth and tilt angles walk through image and project factors onto image voxel
    for( int tilt = 0; tilt < segment_count; tilt++) {      // may want to restrict this to "positive" tilt's
        if( tilt > segment_count/2)
            dz = (tilt-segment_count) * span * crystal_pitch;
        else
            dz = tilt * span * crystal_pitch;
        
        for( int azi = 0; azi < azimuth_count; azi++ ) {
            if(azi > azimuth_count/2)
                dx = (azi - azimuth_count)*span*crystal_pitch;
            else
                dx = azi*span*crystal_pitch;
            
            phi = atan2( dx, detector_distance);           // azimuth angle ( in rads )
            double sinphi = sin(phi);
            double cosphi = cos(phi);
            double tanphi = dx/detector_distance;
            
            temp = detector_distance * detector_distance + dx * dx;
            theta = atan( dz / sqrt(temp));     // tilt (or polar) angle
            double costheta = cos(theta);
            double sintheta = sin(theta);

            int plane_index  = (tilt * azimuth_count + azi) * projection_size;
            int range = round(FWHM)+1;      // FWHM in bin units;
            int left, right;
            int top, bottom;

            for( int k=0; k<zbins; k++) {                       // current voxelsize = 0.8 mm
                z = (k - zbins/2) * voxel_size;                 // z-coordinate of voxel in mm;  zbins/2 = 58;
                for( int j=0; j<ybins; j++) {
                    y = (j-ybins/2) * voxel_size;               // y-coordinate of voxel in mm;  ybins/2 = 22;
                    for(int i=0; i< xbins; i++) {
                        int mask_index = i + j*xbins + k*ybins*xbins;
                        mask_index += (tilt * azimuth_count + azi) * (zbins*ybins*xbins);
                        if(image_mask[mask_index] == 0) continue;
                        x = (i-xbins/2) * voxel_size;           // x-coordinate of voxel in mm;  xbins/2 = 25;
                        
                        r = (-y * sinphi + x * cosphi);         // radial distance (in mm) in projection
                        if( fabs(r) > 20.0) continue;
                        double arg = y*y + x*x - r*r;
                        if(arg < 0.0) {
//                          printf("WARNING: tilt=%2d azi=%2d x = %6.3lf y = %6.3lf z = %6.1lf r %6.3lf and arg %6.3lf\n", tilt, azi, x, y, z, r, arg);
                            continue;
                        }
                        if( y < x*tanphi){
                            s = z * costheta + sqrt(arg) * sintheta;    // axial location (in mm) in projection
                        }
                        else {
                            s = z * costheta - sqrt(arg) * sintheta;
                        }
                        if( fabs(s) > 46.8) continue;
                        r /= bin_width;
                        r += xcenterbins;
                        s /= bin_width;
                        s += zcenterbins;
                        int  ip = floor(r);
                        int  jp = floor(s);
                        
                        // implementing Gaussian blurring here
                        left   = MAX((ip-range), 0);
                        right  = MIN((ip+range), radial_bins);
                        bottom = MAX((jp-range), 0);
                        top    = MIN((jp+range), axial_bins);
        //                        if(ctr < 10) printf("left=%2d  right=%2d  ", left, right);
                        for( int ii=left; ii<right; ii++) {
                            wtx[ii] = gaussian(r-ii,FWHM);
                            if(wtx[ii] < 0.01) {
                                if(ii<ip) left++;
                                if(ii>ip) right = ii;
                            }
                        }
        //                        if(ctr++ < 10) printf("left=%2d  ip=%2d  right=%2d  n=%2d\n", left, ip, right, right - left);

                        float tmp = 0.0;
                        float ttt = 0.0;
                        for( int jj=bottom; jj<top; jj++) {
                            float wty = gaussian(s-jj,FWHM);
                            if(wty < 0.01) continue;
                            for( int ii=left; ii<right; ii++) {
                                int index = jj * radial_bins + ii;
                                tmp += est_planogram[plane_index+index]*wtx[ii]*wty;
                                ttt += wtx[ii]*wty;
                            }
                        }
                        int voxel = i + j*xbins+ k*ybins*xbins;
                        factors[voxel] += tmp;
                        weights[voxel] += ttt;
                    }
                }
            }
        }
    }
    
    // update image (finally!)
    for( int i=0; i< xbins*ybins*zbins; i++) if(weights[i]>0.0) image[i] *= (factors[i]/weights[i]);
    
    return;
}


//************************************** calculate normalization correction **************************
void  calculate_normalization_correction( float *measured_planogram, int planogram_size)
{
    // measured data have been binned into measured_planogram; now we have to correct the various
    // projections for differences in flood phantom thickness according to tilt and azimuth angles
    // First, normalize counts per straight projection (tilt=azimuth=0), then apply angle corrections to tother projections
    
    double sum = 0.0;
    int    index=0;
    for(int j=0; j<2*DETECTOR_ROWS-1; j++) {
        for( int i=0; i< 2*DETECTOR_COLS-1; i++) {
            sum += measured_planogram[index];
            index++;
        }
    }
    double avg_per_pixel = sum/index;
    printf("avg_per_pixel = %10.3lf\n", avg_per_pixel);
    
    double  dx, dz;
    index=0;
    for( int tilt = 0; tilt < segment_count; tilt++) {      // may want to restrict this to "positive" tilt's
        if( tilt > segment_count/2)
            dz = (tilt-segment_count) * span * crystal_pitch;
        else
            dz = tilt * span * crystal_pitch;
        
        for( int azi = 0; azi < azimuth_count; azi++ ) {
            if(azi > azimuth_count/2)
                dx = (azi - azimuth_count)*span*crystal_pitch;
            else
                dx = azi*span*crystal_pitch;
            
            double phi = atan2( dx, detector_distance);     // azimuth angle ( in rads )
            double cosphi = cos(phi);
            
            double temp = detector_distance * detector_distance + dx * dx;
            double theta = atan( dz / sqrt(temp));     // tilt (or polar) angle
            double costheta = cos(theta);
            
            for(int j=0; j<2*DETECTOR_ROWS-1; j++) {
                for( int i=0; i< 2*DETECTOR_COLS-1; i++) {
                    measured_planogram[index] /= (avg_per_pixel/(costheta*cosphi));
                    index++;
                }
            }
        }
    }
    // write result to disk
    write_r4_array( &measured_planogram[0], planogram_size, "normalization_planogram.img");
    
}
//************************************** create measured planogram **************************
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
    int time=0, crystalindex1, crystalindex2, eventblock[EVENTBLOCK_SIZE];
    
    unsigned int rawimageD1[DETECTOR_ROWS * DETECTOR_COLS] = {0};
    unsigned int rawimageD2[DETECTOR_ROWS * DETECTOR_COLS] = {0};
    
    double y1 = (detector_distance * (-0.5));
    double y2 = (detector_distance * 0.5);
    double xcenter = (DETECTOR_COLS - 1) * 0.5;
    double zcenter = (DETECTOR_ROWS - 1) * 0.5;
    double r, phi, s, theta;
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
        
        sino_coord(x1, x2, y1, y2, z1, z2, &r, &phi, &theta, &s);
        r += xcenterbins;
        s += zcenterbins;
        // implementing bilinear interpolation here
        int   i = floor(r);
        int   j = floor(s);
        float fx = r-i;
        float fy = s-j;
        
        int index = j * radial_bins + i;
        index += azimuth_segment * projection_size;
        index += tilt_segment * azimuth_count * projection_size;
        if (index < planogram_size && index >= 0){
            planogram[index]   += (1.0-fx)*(1.0-fy);
            planogram[index+1] +=  fx*(1.0-fy);
            planogram[index+radial_bins]   +=  (1.0-fx)*fy;
            planogram[index+radial_bins+1] +=  fx*fy;
        } else {
            printf("Something's wrong with r= %6.3lf and s= %6.3lf\n", r, round(s));
        }
        event_count++;
        if( (event_count % 1000000)==0 ) printf("%10d events processed. scan time = %7.1lf minutes\n", event_count, (float) time/60000.0);
        
    }
    
    printf("%10d events processed. scan time = %7.1lf minutes\n", event_count, (float) time/60000.0);
    fclose( fp );
    
    char  filename[255];
    //    sprintf(filename,"%s_Det1_raw.img", basename);
    //    write_i4_array( &rawimageD1[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    //    sprintf(filename,"%s_Det2_raw.img", basename);
    //   write_i4_array( &rawimageD2[0], DETECTOR_ROWS*DETECTOR_COLS, filename);
    
    sprintf(filename,"%s_measured_projections.img", basename);
    write_r4_array( &planogram[0], planogram_size, filename);
    
}

//************************************** create_normplanogram **************************
void create_normplanogram( float *normplanogram, const int planogram_size)
{   // we are assigning one count to each line-of-response and bin it into the planogram
    // this can be used in lieu of a true normalization correction
    // Main use of this routine is fo evaluation of binning schemes
    // we can get rid of this routine after initial testing phase is over
    
    char    filename[255];
    double y1 = (detector_distance * (-0.5));
    double y2 = (detector_distance * 0.5);
    double xcenter = (DETECTOR_COLS - 1) * 0.5;
    double zcenter = (DETECTOR_ROWS - 1) * 0.5;     // zcenter = 29.0
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
                printf("row1=%2d row2=%2d col1=%2d \n", row1, row2, col1);
                for (int col2 = 0; col2 < DETECTOR_COLS; col2++) {
                    int azimuth_segment = col2 - col1;  // NEW
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
                    
                    if (index < planogram_size && index >= 0){
                        normplanogram[index]   += (1.0-fx)*(1.0-fy);
                        normplanogram[index+1] +=  fx*(1.0-fy);
                        if(j<116) {
                            normplanogram[index+radial_bins]   +=  (1.0-fx)*fy;
                            normplanogram[index+radial_bins+1] +=  fx*fy;
                        }
                    } else {
                        printf("Something's wrong with r %6.3lf and z %6.3lf\n", r, round(z));
                    }
                }
            }
        }
    }
    printf("norm planogram created.\n");
    sprintf(filename, "test_norm_planogram.img");
    write_r4_array(&normplanogram[0], planogram_size, filename);
}


//************************************** create image mask *************************************************************
void create_image_mask( short int *image_mask, int xbins, int ybins, int zbins, float *normplanogram, int planogram_size)
{
    // note that  xbins, ybins, zbins here a half of image xbins, etc. used in forward_project()
    // voxelsize here is 1.6mm instead of 0.8mm
    
    double xcenter = (DETECTOR_COLS - 1) * 0.5;
    double zcenter = (DETECTOR_ROWS - 1) * 0.5;
    double xcenterbins = xcenter * crystal_pitch /bin_width;
    double zcenterbins = zcenter * crystal_pitch /bin_width;
    double  temp, theta, phi, x,y,z, r,s;
    double dx, dz;
    double mask_voxel_size = voxel_size;
    
    // for all azimuth and tilt angles walk through image and project onto planogram
    for( int tilt = 0; tilt < segment_count; tilt++) {      // may want to restrict this to "positive" tilt's
        if( tilt > segment_count/2)
            dz = (tilt-segment_count) * span * crystal_pitch;
        else
            dz = tilt * span * crystal_pitch;
        
        for( int azi = 0; azi < azimuth_count; azi++ ) {
            if(azi > azimuth_count/2)
                dx = (azi - azimuth_count)*span*crystal_pitch;
            else
                dx = azi*span*crystal_pitch;
            
            phi = atan2( dx, detector_distance);           // azimuth angle ( in rads )
            double sinphi = sin(phi);
            double cosphi = cos(phi);
            double tanphi = dx/detector_distance;
            
            temp = detector_distance * detector_distance + dx * dx;
            theta = atan( dz / sqrt(temp));     // tilt (or polar) angle
            double costheta = cos(theta);
            double sintheta = sin(theta);
            
            // xbins, ybins, zbins here a half of image xbins, etc.
            for( int k=0; k<zbins; k++) {                       // this voxelsize = 1.6 mm
                z = (k - zbins/2) * mask_voxel_size;                 // z-coordinate of voxel in mm
                for( int j=0; j<ybins; j++) {
                    y = (j-ybins/2) * mask_voxel_size;               // y-coordinate of voxel in mm
                    for( int i=0; i< xbins; i++) {
                        x = (i-xbins/2) * mask_voxel_size;           // x-coordinate of voxel in mm
                        r = (-y * sinphi + x * cosphi);         // radial distance (in mm) in projection
                        double arg = y*y + x*x - r*r;
                        if(arg < 0.0) {
                            if(arg > -5.0E-13)
                                arg = 0.0;
                            else {
                                printf("WARNING: tilt=%2d azi=%2d sinphi=%8.4lf cosphi=%8.4lf x= %10.4lf y= %10.4lf z= %7.2lf r=%10.6lf and arg= %13.9lf\n", tilt, azi, sinphi, cosphi, x, y, z, r, arg);
                                continue;
                            }
                        }
                        if( y < x*tanphi){
                            s = z * costheta + sqrt(arg) * sintheta;    // axial location (in mm) in projection
                        }
                        else {
                            s = z * costheta - sqrt(arg) * sintheta;
                        }
                        r /= bin_width;
                        r += xcenterbins;
                        s /= bin_width;
                        s += zcenterbins;
                        
                        int  ip = round(r);
                        int  jp = round(s);
                        if( ip<0 || ip > 50) continue;      // out of projection
                        if( jp<0 || jp > 116) continue;      // out of projection
                        int index = jp * radial_bins + ip;      //projection_size = radial_bins * axial_bins;
                        
                        index += (tilt * azimuth_count + azi) * projection_size;
                        if(normplanogram[index]>0.01) {
                            int maskvoxel = i + j*xbins+ k*ybins*xbins;
                            maskvoxel += (tilt * azimuth_count + azi) * (zbins*ybins*xbins);
                            //                            printf("        maskvoxel2 = %8d\n", maskvoxel);
                            image_mask[ maskvoxel ]++;
                        }
                    }
                }
            }
        }
    }
    int image_mask_size = segment_count*azimuth_count * (xbins*ybins*zbins);
    write_i2_array( image_mask, image_mask_size, "image_mask.raw");
}


//************************************** forward_projection of estimated image **************************
void forward_project( float *image, short int *image_mask, int xbins, int ybins, int zbins, float *est_planogram, int planogram_size)
{
    double xcenter = (DETECTOR_COLS - 1) * 0.5;
    double zcenter = (DETECTOR_ROWS - 1) * 0.5;
    double xcenterbins = xcenter * crystal_pitch /bin_width;        // = 25
    double zcenterbins = zcenter * crystal_pitch /bin_width;        // = 58
    double temp, theta, phi, x,y,z, r,s;
    double dx, dz, FWHM = 1.5;      // FWHM in bin_width units;  1.5 corresponds to 1.2 mm
    int    ctr = 0;
    float  wtx[radial_bins];
    
    // for all azimuth and tilt angles walk through image and project onto planogram
    for( int tilt = 0; tilt < segment_count; tilt++) {      // may want to restrict this to "positive" tilt's
        if( tilt > segment_count/2)
            dz = (tilt-segment_count) * span * crystal_pitch;
        else
            dz = tilt * span * crystal_pitch;
        
        for( int azi = 0; azi < azimuth_count; azi++ ) {
            if(azi > azimuth_count/2)
                dx = (azi - azimuth_count)*span*crystal_pitch;
            else
                dx = azi*span*crystal_pitch;
            
            phi = atan2( dx, detector_distance);           // azimuth angle ( in rads )
            double sinphi = sin(phi);
            double cosphi = cos(phi);
            double tanphi = dx / detector_distance;
            
            temp = detector_distance * detector_distance + dx * dx;
            theta = atan( dz / sqrt(temp));     // tilt (or polar) angle
            double costheta = cos(theta);
            double sintheta = sin(theta);
            
            int plane_index  = (tilt * azimuth_count + azi) * projection_size;
            int range = round(FWHM)+1;      // FWHM in bin units;
            int left, right;
            int top, bottom;
            
            // loop over image space
            for( int k=0; k<zbins; k++) {                       // current voxelsize = 0.8 mm
                z = (k - zbins/2) * voxel_size;                 // z-coordinate of voxel in mm;  zbins/2 = 58;
                for( int j=0; j<ybins; j++) {
                    y = (j-ybins/2) * voxel_size;               // y-coordinate of voxel in mm;  ybins/2 = 22;
                    for(int i=0; i< xbins; i++) {
                        int mask_index = i + j*xbins + k*ybins*xbins;
                        mask_index += (tilt * azimuth_count + azi) * (zbins*ybins*xbins);
                        if(image_mask[mask_index] == 0) continue;
                        x = (i-xbins/2) * voxel_size;           // x-coordinate of voxel in mm;  xbins/2 = 25;
                        
                        r = (-y * sinphi + x * cosphi);         // radial distance (in mm) in projection
                        if( fabs(r) > 20.0) {
                            printf("ERROR: tilt=%2d azi=%2d phi = %7.3lf, x = %6.3lf y = %6.3lf  z = %6.1lf  r = %6.3lf\n", tilt, azi, phi, x, y, z, r);
                            if(ctr++ > 10) exit(-1);
                            continue;
                        }
                        double arg = y*y + x*x - r*r;
                        if(arg < 0.0) {
                            if(arg > -5.0E-13) {
                                arg = 0.0;
                            }
                            else {
                                printf("WARNING: tilt=%2d azi=%2d x = %6.3lf y = %6.3lf  z = %6.1lf  r %6.3lf and arg %6.3lf\n", tilt, azi, x, y, z, r, arg);
                                continue;
                            }
                        }
                        if( y < x*tanphi){
                            s = z * costheta + sqrt(arg) * sintheta;    // axial location (in mm) in projection
                        }
                        else {
                            s = z * costheta - sqrt(arg) * sintheta;
                        }
                        if( fabs(s) > 46.8) continue;
                        r /= bin_width;
                        r += xcenterbins;
                        s /= bin_width;
                        s += zcenterbins;
                        int  ip = floor(r);
                        int  jp = floor(s);

                        // implementing Gaussian blurring here
                        left   = MAX((ip-range), 0);
                        right  = MIN(ip+range,radial_bins);
                        bottom = MAX(jp-range, 0);
                        top    = MIN(jp+range, axial_bins);
//                        if(ctr < 10) printf("left=%2d  right=%2d  ", left, right);
                        for( int ii=left; ii<right; ii++) {
                            wtx[ii] = gaussian(r-ii,FWHM);
                            if(wtx[ii] < 0.01) {
                                if(ii<ip) left++;
                                if(ii>ip) right = ii;
                            }
                        }
//                        if(ctr++ < 10) printf("left=%2d  ip=%2d  right=%2d  n=%2d\n", left, ip, right, right - left);

                        int voxel = i + j*xbins+ k*ybins*xbins;
                        float counts  = image[voxel];
                        for( int jj=bottom; jj<top; jj++) {
                            float wty = gaussian(s-jj,FWHM);
                            if(wty < 0.01) continue;
                            for( int ii=left; ii<right; ii++) {
                                int index = jj * radial_bins + ii;
                                est_planogram[plane_index+index] += counts*wtx[ii]*wty;
                            }
                        }
                    }
                }
            }
            
        }
    }
}

//************************************** Gaussian *******************************************
float gaussian( double x, double FWHM)
{
    float result = exp(-2.772589*(x*x/(FWHM*FWHM)))/(1.06447*FWHM);
    return result;
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
void sino_coord( double x1, double x2, double y1, double y2, double z1, double z2, double* r, double* phi, double* theta, double *s)
{
    *phi = atan2(x2 - x1, y2 - y1);           // azimuth angle ( in rads )
    *r = (-y2 * sin(*phi) + x2 * cos(*phi));  // radial distance (in mm)
    //    *phi *= (180.0 / M_PI) / angular_bin_width;   // converts to angular bin number  (not needed currently)
    
    double temp = (y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1);
    *theta = atan((z2 - z1) / sqrt(temp));     // tilt (or polar) angle
    
    double arg = y1 * y1 + x1 * x1 - (*r)*(*r);
    *s =sqrt(arg) * sin(*theta) + z1 * cos(*theta);  // axial location (in mm)
    *s /= bin_width;
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
void parse_command_line(int argc, const char *argv[], char inputfile[], char basename[], char path[])
{
    // command line must specify list file name first and it must be followed
    //    by the other parameters
    // e.g.  ./listproc LYSOtest.bin -a400 -b650
    
    if(argc <= 1) {
        terminate_with_error("Not enough arguments to find filename!");
    }
    
    if( sscanf( &argv[1][0], "%[^\t\n]", &inputfile[0]) <= 0 ) {
        terminate_with_error("No luck parsing command line for filename!");
    }
    
    printf("Filename = %s\n", inputfile);
    // list file or projections?
    // basename = list file name without extension
    char *filename;
    filename = strrchr(inputfile, '/');
    if(filename != NULL) {
        size_t len = strlen(inputfile) - strlen(filename);
        strncpy( path, inputfile, len);
    }
    else {
        sprintf(path, "");
    }
    
    char * extension = strstr(inputfile,".bin");
    if(extension == NULL) {  // projection file?
        extension = strstr(inputfile,"_measured_projections.img");
        if(extension == NULL) {
            terminate_with_error("No luck parsing command line for proper extension!");
        }
        else {
            size_t len = strlen(inputfile) - strlen(extension);
            strncpy( basename, inputfile, len);
        }
    }
    else {
        // len is the length of the basename
        size_t len = strlen(inputfile) - strlen(extension);
        strncpy( basename, inputfile, len);
        is_listfile = 1;
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
                case ITERATIONS_FLAG:
                    iterations = atoi(&argv[i][2]);
                    printf("doing %d iterations!\n", iterations);
                    break;
                case NORMALIZATION_FLAG:
                    process_normalization = 1;
                    printf("creating normalization planogram.\n");
                    break;
                case TEST_FLAG:
                    testing = 1;
                    printf("testing testing testing\n");
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
void parse_from_stdin( char listfile[], char basename[], char path[] )
{
    // get arguments interactively
    char buf[255];
    printf( "Enter listfile path:\n" );
    fgets( buf, 255, stdin );
    
    if ( sscanf( &buf[0], "%[^\t\n]", &listfile[0]) > 0 ) {
        
        char *filename = strrchr(listfile, '/');
        size_t len = strlen(listfile) - strlen(filename);
        strncpy( path, listfile, len);
        
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
void write_i2_array(short array[], int nelements, char filename[])
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

//************************************** read_r4_array ********************
void read_r4_array(float array[], int nelements, char filename[])
{
    FILE *fp = fopen( filename,"rb");
    if ( !fp ) {
        char message[255];
        sprintf(message,"Error opening file %s\n", filename);
        terminate_with_error(message);
    }
    
    size_t n = fread( &array[0], 4, nelements, fp);
    if( n != nelements ) {
        char message[255];
        sprintf(message,"Write error! Only %zd words written instead of %d, file = %s\n", n, nelements, filename);
        terminate_with_error(message);
    }
    
    printf ("File %s was read from disk.\n", filename);
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
