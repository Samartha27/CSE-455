#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
    float x1 = ((x - floor(x)) >0.5) ?ceil(x):floor(x+0.5);
    float y1 = ((y - floor(y)) >0.5) ?ceil(y):floor(y+0.5);

    return get_pixel(im,x1,y1,c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    image im_new = make_image(w,h,im.c);
    //printf("I am here");

    float w_ratio = (float)im.w/w; // a_width(w-0.5) + b_width = (im.w-0.5)
    float h_ratio = (float)im.h/h;// old y /new y = old height/new height

    for (int c = 0;c<im.c;c++){
        for (int y = 0;y<h;y++){
            for (int x = 0;x<w;x++){
                float x_old = (x * w_ratio) + 0.5*(w_ratio -1);
                float y_old = (y * h_ratio) + 0.5*(h_ratio -1); 

                float old_pixel_value = nn_interpolate(im,x_old,y_old,c);
                set_pixel(im_new,x,y,c,old_pixel_value);
            }
        }
    }
    return im_new;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO In the old image get the coordinates of the pixel which falls inbetween give set of 4 pixels
    float x1_y1_v = get_pixel(im,floorf(x),floorf(y),c);
    float x2_y1_v = get_pixel(im,ceilf(x),floorf(y),c);
    float x1_y2_v = get_pixel(im,floorf(x),ceilf(y),c);
    float x2_y2_v = get_pixel(im,ceilf(x),ceilf(y),c);

    float d1 = x - floorf(x);
    float d2 = ceilf(x) - x;
    float d3 = y - floorf(y);
    float d4 = ceilf(y) - y;

    float q1 = x1_y1_v*d2 + x2_y1_v*d1;
    float q2 = x1_y2_v*d2 + x2_y2_v*d1;
    float q = q1*d4 + d3*q2;

    return q;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    //printf("I am here");
    image im_new =  make_image(w,h,im.c);

    float w_ratio = (float) im.w/w;
    float h_ratio = (float) im.h/h;
    

    for (int c = 0;c<im.c;c++){
        for (int y = 0;y<h;y++){
            for (int x = 0;x<w;x++){
                //printf("I am inside loop");
                float x_old = (x * w_ratio) +  0.5*(w_ratio -1);
                float y_old = (y * h_ratio) +  0.5*(h_ratio -1);  // old y /new y = old height/new height
                
                float old_pixel_value = bilinear_interpolate(im,x_old,y_old,c);
                set_pixel(im_new,x,y,c,old_pixel_value);
            }
        }
    }
    //printf("I am outside loop");

    return im_new;
}

