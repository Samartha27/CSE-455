#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    float sum = 0;
        for (int y = 0;y<im.h;y++){
            for (int x = 0;x<im.w;x++){
                sum +=get_pixel(im, x, y, 0);
            }}
    if (sum ==0){
        return;
    }
    float value = 0.0;
    
        for (int y = 0;y<im.h;y++){
            for (int x = 0;x<im.w;x++){
                value = get_pixel(im, x, y, 0)/sum;
                set_pixel(im,x,y,0,value);
            }}


}

image make_box_filter(int w)
{
    // TODO
    image box_filter = make_image(w,w,1);
    
        for (int y = 0;y<box_filter.w;y++){
            for (int x = 0;x<box_filter.w;x++){
                set_pixel(box_filter,x,y,0,1.0/(w*w));
            }}
    return box_filter;
}


image convolve_image(image im, image filter, int preserve)
{
    // TODO
    assert(filter.c == im.c || filter.c ==1);
    image im_new = make_image(im.w,im.h,im.c);
    
    
        for (int c = 0;c<im_new.c;c++){
            for (int y = 0;y<im_new.h;y++){
                for (int x = 0;x<im_new.w;x++){
                    float new_value =0.0;
                    for (int j = 0;j<filter.h;j++){
                        for (int i = 0;i<filter.w;i++){
                            int x_d = filter.w/2;
                            int y_d = filter.h/2;
                            new_value += get_pixel(im,x+i-x_d,y+j-y_d,c)*get_pixel(filter,i,j,filter.c==1?0:c);
                        }}
                        set_pixel(im_new,x,y,c,new_value);
                }}}

    if (preserve!=1){
        image im_pnew = make_image(im_new.w,im_new.h,1);
    
        for (int y = 0;y<im_new.h;y++){
            for (int x = 0;x<im_new.w;x++){
                float new_value = 0.0;
                for(int pc=0;pc<im_new.c;pc++){
                    new_value += get_pixel(im_new,x,y,pc);
                }
                set_pixel(im_pnew,x,y,0,new_value);
            }}
            return im_pnew;
    }

    return im_new;
}

image make_highpass_filter()
{
    // TODO
    image high_pass = make_box_filter(3);
    high_pass.data[0] = 0;
    high_pass.data[1] = -1;
    high_pass.data[2] = 0;
    high_pass.data[3] = -1;
    high_pass.data[4] = 4;
    high_pass.data[5] = -1;
    high_pass.data[6] = 0;
    high_pass.data[7] = -1;
    high_pass.data[8] = 0;

    return high_pass;
}

image make_sharpen_filter()
{
    // TODO
    image sharpen = make_box_filter(3);
    float b[9] = {0,-1,0,-1,5,-1,0,-1,0};

    memcpy(sharpen.data,b,sizeof(b));

    return sharpen;
    
}

image make_emboss_filter()
{
    // TODO
    image emboss = make_image(3,3,1);
    float c[9] = {-2,-1,0,-1,1,1,0,1,2};

    memcpy(emboss.data,c,sizeof(c));

    return emboss;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: We use preserve for box filter,emboss filter and sharpen filter because we want to retain all the channels from the  original image
// with only certain specific transformations in order to highlight or lighten certain characteristics.
// We do not use preserve for high pass filter since we only require a certain subset of the original image 
//(the edges from the original image// wherein just one channel is sufficient to show the edges 
//and red , green and blue channels are not necessary.

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Yes. When a high pass filter is used,
// we will want to do postprocessing with the help of clamp to limit the overflow of pixel values in the image
// so that the edge of the dog is visible clearly.

image make_gaussian_filter(float sigma)
{
    // TODO
    int w = (int) roundf(6*sigma)+1;
    
    image gaussian_filter = make_box_filter(w);
    int offset = (int) w/2;
        for (int y = 0;y<w;y++){
            for (int x = 0;x<w;x++){
                float g = (1/(TWOPI*sigma*sigma))*expf(-(pow(x-offset,2)+pow(y-offset,2))/(2*pow(sigma,2)));
                
                set_pixel(gaussian_filter,x,y,0,g);
            }}
    l1_normalize(gaussian_filter);
    return gaussian_filter;
}

image add_image(image a, image b)
{
    // TODO
    if (a.w != b.w || a.h != b.h){
        return make_image(1,1,1);
    }
    image im_new = make_image(a.w,a.h,a.c);
    for (int c = 0;c<im_new.c;c++){
            for (int y = 0;y<im_new.h;y++){
                for (int x = 0;x<im_new.w;x++){
                    float value = get_pixel(a,x,y,c)+get_pixel(b,x,y,c);
                    set_pixel(im_new,x,y,c,value);
                }}}

    return im_new;
}

image sub_image(image a, image b)
{
    // TODO
    if (a.w != b.w || a.h != b.h ){
        return make_image(1,1,1);
        }
    image im_new = make_image(a.w,a.h,a.c);
    for(int x=0;x<a.w*a.h*a.c;x++){
        im_new.data[x] = a.data[x] - b.data[x];
    }
    return im_new;
}

image make_gx_filter()
{
    // TODO
    image gx_filter = make_box_filter(3);
    float d[9] = {-1,0,1,-2,0,2,-1,0,1};

    memcpy(gx_filter.data,d,sizeof(d));

    return gx_filter;
}

image make_gy_filter()
{
    // TODO
    image gy_filter = make_box_filter(3);
    float e[9] = {-1,-2,-1,0,0,0,1,2,1};

    memcpy(gy_filter.data,e,sizeof(e));

    return gy_filter;
}

void feature_normalize(image im)
{
    // TODO
    float min = im.data[0];
    float max = im.data[0];

    for(int i =0;i<im.w*im.h*im.c;i++){
        if (min > im.data[i]){min = im.data[i];}
        if (max < im.data[i]){max = im.data[i];}
    }
    float range = max- min;
    if (range==0){
        for(int i =0;i<im.w*im.h*im.c;i++){
            im.data[i]=0;  }}
    else {
        for(int i =0;i<im.w*im.h*im.c;i++){
            im.data[i]=(im.data[i]-min)/range; 
        }}
}
    

image *sobel_image(image im)
{
    // TODO
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gx_transform = convolve_image(im,gx_filter,0);
    image gy_transform = convolve_image(im,gy_filter,0);

    image grad_mag = make_image(im.w,im.h,1);
    image grad_dir = make_image(im.w,im.h,1);

    for(int i=0;i<im.w*im.h;i++){
        grad_mag.data[i] = sqrtf(pow(gx_transform.data[i],2) +pow(gy_transform.data[i],2));
        grad_dir.data[i] = atan2f(gy_transform.data[i],gx_transform.data[i]);
    }

    image *sobel;
    sobel = calloc(2, sizeof(image));
    *(sobel+0)=grad_mag;
    *(sobel+1)=grad_dir;
    //sobel[0]= grad_mag;
    //sobel[1]= grad_dir;
    return sobel;
}

image colorize_sobel(image im)
{
    // TODO
    image hsv_im = make_image(im.w,im.h,3);
    image *sobel_im = calloc(2, sizeof(image));
    sobel_im = sobel_image(im);
    image sob_mag = sobel_im[0];
    image sob_dir = sobel_im[1];
    feature_normalize(sob_mag);
    feature_normalize(sob_dir);
        for (int y = 0;y<im.h;y++){
            for (int x = 0;x<im.w;x++){
                set_pixel(hsv_im,x,y,0,get_pixel(sob_dir,x,y,0));
                set_pixel(hsv_im,x,y,1,get_pixel(sob_mag,x,y,0));
                set_pixel(hsv_im,x,y,2,get_pixel(sob_mag,x,y,0));
            }}
    hsv_to_rgb(hsv_im);
    image f = make_gaussian_filter(2);
    hsv_im = convolve_image(hsv_im,f,1);

    return hsv_im;
}
