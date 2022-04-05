#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <unistd.h> 

float get_pixel_zero_pad(image im, int x, int y, int c)
{
	//printf("%d, %d \n",x, y, c);
	if(x < 0){
		return 0.0;
	}
	if(y < 0){
		return 0.0;
	}
	if(x >= im.w){
		return 0.0;
	}
	if(y >= im.h){
		return 0.0;
	}
    long index = c*im.h*im.w + y*im.w + x;
    float returnVal = im.data[index];
    return returnVal;
}

// Draws a line on an image with color corresponding to the direction of line
// image im: image to draw line on
// float x, y: starting point of line
// float dx, dy: vector corresponding to line angle and magnitude
void draw_line(image im, float x, float y, float dx, float dy)
{
    assert(im.c == 3);
    float angle = 6*(atan2(dy, dx) / TWOPI + .5);
    int index = floor(angle);
    float f = angle - index;
    float r, g, b;
    if(index == 0){
        r = 1; g = f; b = 0;
    } else if(index == 1){
        r = 1-f; g = 1; b = 0;
    } else if(index == 2){
        r = 0; g = 1; b = f;
    } else if(index == 3){
        r = 0; g = 1-f; b = 1;
    } else if(index == 4){
        r = f; g = 0; b = 1;
    } else {
        r = 1; g = 0; b = 1-f;
    }
    float i;
    float d = sqrt(dx*dx + dy*dy);
    for(i = 0; i < d; i += 1){
        int xi = x + dx*i/d;
        int yi = y + dy*i/d;
        set_pixel(im, xi, yi, 0, r);
        set_pixel(im, xi, yi, 1, g);
        set_pixel(im, xi, yi, 2, b);
    }
}

// Make an integral image or summed area table from an image
// image im: image to process
// returns: image I such that I[x,y] = sum{i<=x, j<=y}(im[i,j])

image make_integral_image(image im)
{
    float sum;
    int x,y,z;
    image integ = make_image(im.w, im.h, im.c);
    // TODO: fill in the integral image
    for ( z = 0;z<im.c;z++){
        for ( y = 0;y<im.h;y++){
            for ( x = 0;x<im.w;x++){
                sum = 0.0;
                sum += get_pixel(im,x,y,z);
                if (x> 0){
                    sum += get_pixel(integ,x-1,y,z);
                }
                if ( y>0){
                    sum += get_pixel(integ,x,y-1,z);
                }
                if (x>0 && y>0){
                    sum -= get_pixel(integ,x-1,y-1,z);
                }
                set_pixel(integ,x,y,z,sum);
            }
        }
    }
    //save_image(integ,"integ");
    return integ;
}


// Apply a box filter to an image using an integral image for speed
// image im: image to smooth
// int s: window size for box filter
// returns: smoothed image
image box_filter_image(image im, int s)
{
    int i,j,k,i_A,i_D,j_A,j_D,window_i,window_j;
    float v1,v2,v3,v4,value,final_value;
    image integ = make_integral_image(im);
    image S = make_image(im.w, im.h, im.c);
    // TODO: fill in S using the integral image.
    int w =  (int) (s/2);
    for ( k = 0;k<im.c;k++){
        for ( j = 0;j<im.h;j++){
            for ( i = 0;i<im.w;i++){
                v1=0.0;v2=0.0;v3=0.0;v4=0.0;

                i_A = i-w-1;
                j_A = j-w-1;
                i_D = i+w;
                j_D = j+w;

                if(i_A<0) i_A=-1;
                if(j_A<0) j_A=-1;
                if(i_D>=im.w) i_D=im.w-1;
                if(j_D>=im.h) j_D=im.h-1;
                
                v4 = get_pixel(integ, i_D, j_D, k);
                if(i_A>=0) v3 = get_pixel(integ, i_A, j_D, k);
                if(j_A>=0) v2 = get_pixel(integ, i_D, j_A, k);
                if((i_A>=0)&&(j_A>=0)) v1 = get_pixel(integ, i_A, j_A, k);

                window_i = i_D - i_A;
                window_j = j_D - j_A;

                value = v1+v4-v3-v2;
                final_value = value/window_i/window_j;

                set_pixel(S,i,j,k,final_value);
            }
        }
    }
    //save_image(S,"box_filter_integral");
    return S;
}

// Calculate the time-structure matrix of an image pair.
// image im: the input image.
// image prev: the previous image in sequence.
// int s: window size for smoothing.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          3rd channel is IxIy, 4th channel is IxIt, 5th channel is IyIt.
image time_structure_matrix(image im, image prev, int s)
{
    int converted = 0;
    if(im.c == 3){
        converted = 1;
        im = rgb_to_grayscale(im);
        prev = rgb_to_grayscale(prev);
    }

    float Ix, Iy, It, Ix2, Iy2;

    image filter_x = make_gx_filter();
    image filter_y = make_gy_filter();

    image Gx = convolve_image(im, filter_x, 0);
    image Gy = convolve_image(im, filter_y, 0);

    //feature_normalize(Gx);
    //save_image(Gx, "Gx");

    image S = make_image(im.w, im.h, 5);

    for (int j = 0;j<im.h;j++){
        for (int i = 0;i<im.w;i++){            
            Ix = get_pixel(Gx, i, j, 0);
            Ix2 = pow(Ix, 2);
            Iy = get_pixel(Gy, i, j, 0);
            Iy2 = pow(Iy, 2);

            It = get_pixel(im, i, j, 0) - get_pixel(prev, i, j, 0);

            set_pixel(S, i, j, 0, Ix2);
            set_pixel(S, i, j, 1, Iy2);
            set_pixel(S, i, j, 2, Ix*Iy);
            set_pixel(S, i, j, 3, Ix*It);
            set_pixel(S, i, j, 4, Iy*It);
        }
    }

    S = box_filter_image(S,s);

    // TODO: calculate gradients, structure components, and smooth them

    if(converted){
        free_image(im); free_image(prev);
    }
    return S;
}

// Calculate the velocity given a structure image
// image S: time-structure image
// int stride: only calculate subset of pixels for speed
image velocity_image(image S, int stride)
{
    image v = make_image(S.w/stride, S.h/stride, 3);
    int i, j;
    matrix M = make_matrix(2,2);
    matrix b = make_matrix(2,1);
    for(j = (stride-1)/2; j < S.h; j += stride){
        for(i = (stride-1)/2; i < S.w; i += stride){
            
            float Ixx = S.data[i + S.w*j + 0*S.w*S.h];
            float Iyy = S.data[i + S.w*j + 1*S.w*S.h];
            float Ixy = S.data[i + S.w*j + 2*S.w*S.h];
            float Ixt = S.data[i + S.w*j + 3*S.w*S.h];
            float Iyt = S.data[i + S.w*j + 4*S.w*S.h];
            
            // TODO: calculate vx and vy using the flow equation

            M.data[0][0] = Ixx;
            M.data[0][1] = Ixy;
            M.data[1][0] = Ixy;
            M.data[1][1] = Iyy;
            //print_matrix(M);
            //printf("I am not here\n");
            matrix M_inv = matrix_invert(M);
            
            b.data[0][0] = -Ixt;
            b.data[1][0] = -Iyt;
            //print_matrix(M);

            if(M_inv.data){
                matrix V = matrix_mult_matrix(M_inv, b);

                set_pixel(v, i/stride, j/stride, 0, V.data[0][0]);
                set_pixel(v, i/stride, j/stride, 1, V.data[1][0]);
                free_matrix(V);
            }
            else{
                //printf("%d, %d\n", i/stride, j/stride);
                set_pixel(v, i/stride, j/stride, 0, 0.0);
                set_pixel(v, i/stride, j/stride, 1, 0.0);
            }
            free_matrix(M_inv);
        }
    }
    free_matrix(M);
    free_matrix(b);

    return v;
}

// Draw lines on an image given the velocity
// image im: image to draw on
// image v: velocity of each pixel
// float scale: scalar to multiply velocity by for drawing
void draw_flow(image im, image v, float scale)
{
    int stride = im.w / v.w;
    int i,j;
    for (j = (stride-1)/2; j < im.h; j += stride) {
        for (i = (stride-1)/2; i < im.w; i += stride) {
            float dx = scale*get_pixel(v, i/stride, j/stride, 0);
            float dy = scale*get_pixel(v, i/stride, j/stride, 1);
            if(fabs(dx) > im.w) dx = 0;
            if(fabs(dy) > im.h) dy = 0;
            draw_line(im, i, j, dx, dy);
        }
    }
}


// Constrain the absolute value of each image pixel
// image im: image to constrain
// float v: each pixel will be in range [-v, v]
void constrain_image(image im, float v)
{
    int i;
    for(i = 0; i < im.w*im.h*im.c; ++i){
        if (im.data[i] < -v) im.data[i] = -v;
        if (im.data[i] >  v) im.data[i] =  v;
    }
}

// Calculate the optical flow between two images
// image im: current image
// image prev: previous image
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// returns: velocity matrix
image optical_flow_images(image im, image prev, int smooth, int stride)
{
    image S = time_structure_matrix(im, prev, smooth);   
    image v = velocity_image(S, stride);
    //constrain_image(v, 6);
    image vs = smooth_image(v, 2);
    free_image(v);
    free_image(S);
    return vs;
}

// Run optical flow demo on webcam
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// int div: downsampling factor for images from webcam
void optical_flow_webcam(int smooth, int stride, int div)
{
#ifdef OPENCV
    void * cap;
    cap = open_video_stream(0, 0, 1280, 720, 30);
    image prev = get_image_from_stream(cap);
    image prev_c = nn_resize(prev, prev.w/div, prev.h/div);
    image im = get_image_from_stream(cap);
    image im_c = nn_resize(im, im.w/div, im.h/div);
    while(im.data){
        image copy = copy_image(im);
        image v = optical_flow_images(im_c, prev_c, smooth, stride);
        draw_flow(copy, v, smooth*div);
        int key = show_image(copy, "flow", 5);
        free_image(v);
        free_image(copy);
        free_image(prev);
        free_image(prev_c);
        prev = im;
        prev_c = im_c;
        if(key != -1) {
            key = key % 256;
            printf("%d\n", key);
            if (key == 27) break;
        }
        im = get_image_from_stream(cap);
        im_c = nn_resize(im, im.w/div, im.h/div);
    }
#else
    fprintf(stderr, "Must compile with OpenCV\n");
#endif
}