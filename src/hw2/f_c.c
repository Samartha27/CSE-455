#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    float l1_norm = 0.0;
    float pix_val;

    for(int i = 0; i<im.h; i++){
        for(int j = 0; j<im.w; j++){
            l1_norm = l1_norm + get_pixel(im, j, i, 0);
        }
    }

    printf("l1_norm - %f\n", l1_norm);

    for(int i = 0; i<im.h; i++){
        for(int j = 0; j<im.w; j++){
            pix_val = get_pixel(im, j, i, 0)/l1_norm;
            set_pixel(im, j, i, 0, pix_val);
            //printf("%f\n", pix_val);
        }
    }

    //printf("im.w - %d, im.h - %d\n", im.w, im.h);
    //float pix_val = 1.0/im.w/im.h;
    //for(int i = 0; i<im.h; i++){
        //for(int j = 0; j<im.w; j++){
            //set_pixel(im, j, i, 0, pix_val);
            //printf("pix_val - %f\n", pix_val);
        //}
    //}
}

image make_box_filter(int w)
{
    image im;
    im = make_image(w, w, 1);
    for(int i = 0; i<im.h; i++){
        for(int j = 0; j<im.w; j++){
            set_pixel(im, j, i, 0, 1);
        }
    }
    l1_normalize(im);
    
    return im;
}

void convolve_layer(image new_im, image im, image filter, int im_c, int filter_c){
    float im_val, fil_val, replace_val;
    int im_x, im_y;

    for(int j = 0; j<im.h; j++){
        for(int k = 0; k<im.w; k++){
            replace_val = 0.0;
            for(int y = 0; y<filter.h; y++){
                for(int x = 0; x< filter.w; x++){
                    im_x = k + x - filter.w/2;
                    im_y = j + y - filter.h/2;
                    
                    im_val = get_pixel(im, im_x, im_y, im_c);
                    fil_val = get_pixel(filter, x, y, filter_c);
                    replace_val = replace_val + im_val*fil_val;

                    //printf("x - %d, y - %d, im_x - %d, im_y - %d, im_val - %f, fil_val - %f, replace_val - %f\n", x, y, im_x, im_y, im_val, fil_val, replace_val);
                }
            }
            set_pixel(new_im, k, j, im_c, replace_val);
            //printf("%f, %f, %f\n", im_val, fil_val, replace_val);
            //if (k == 2) exit(0);
        }
    }
}

image flatten_image(image im){
    image flat_image = make_image(im.w, im.h, 1);
    float val;
    for(int i = 0; i<im.w; i++){
        for(int j = 0; j<im.h; j++){
            val = 0.0;
            for(int k = 0; k<im.c; k++){
                val = val + get_pixel(im, i, j, k);
            }
            //printf("val - %f\n", val);
            set_pixel(flat_image, i, j, 0, val);
        }
    }
    return flat_image;
}

image convolve_image(image im, image filter, int preserve)
{
    assert((im.c == filter.c)||(filter.c == 1));
    image new_im = make_image(im.w, im.h, im.c);
    if(filter.c == 1){
        for(int z = 0; z < im.c; z++){
            convolve_layer(new_im, im, filter, z, 0);
        }
    }
    else{
        for(int z = 0; z < filter.c; z++){
            convolve_layer(new_im, im, filter, z, z);
        }
    }
    if(preserve == 0){
        return flatten_image(new_im);
    }
    return new_im;
}

image make_highpass_filter()
{
    image im = make_image(3,3,1);
    im.data[0] = 0;
    im.data[1] = -1;
    im.data[2] = 0;
    im.data[3] = -1;
    im.data[4] = 4;
    im.data[5] = -1;
    im.data[6] = 0;
    im.data[7] = -1;
    im.data[8] = 0;
    return im;
}

image make_sharpen_filter()
{
    image im = make_image(3,3,1);
    im.data[0] = 0;
    im.data[1] = -1;
    im.data[2] = 0;
    im.data[3] = -1;
    im.data[4] = 5;
    im.data[5] = -1;
    im.data[6] = 0;
    im.data[7] = -1;
    im.data[8] = 0;
    return im;
}

image make_emboss_filter()
{
    image im = make_image(3,3,1);
    im.data[0] = -2;
    im.data[1] = -1;
    im.data[2] = 0;
    im.data[3] = -1;
    im.data[4] = 1;
    im.data[5] = 1;
    im.data[6] = 0;
    im.data[7] = 1;
    im.data[8] = 2;
    return im;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    image im;
    im = make_image(6*sigma + 1, 6*sigma + 1, 1);

    float ret_val, exp_val;
    int x, y;

    for(int im_x = 0; im_x < im.w; im_x++){
        for(int im_y = 0; im_y < im.h; im_y++){
            x = im_x - im.w/2;
            y = im_y - im.h/2;
            ret_val = 1.0/TWOPI/pow(sigma,2);
            exp_val = -(pow(x,2) + pow(y,2))/2.0/pow(sigma,2);
            ret_val = ret_val*exp(exp_val);
            //printf("ret_val - %f, x - %d, y - %d, im_x - %d, im_y - %d\n", ret_val, x, y, im_x, im_y);
            set_pixel(im, im_x, im_y, 0, ret_val);
        }
    }
    l1_normalize(im);
    return im;
}

image add_image(image a, image b)
{
    float val1, val2;
    if((a.h != b.h)&&(a.w != b.w)){
        return make_image(1,1,1);
    }
    image new_im;
    new_im = make_image(a.w, a.h, a.c);
    for(int i = 0; i<a.c; i++){
        for(int j = 0; j<a.h; j++){
            for(int k = 0; k<a.w; k++){
                val1 = get_pixel(a, k, j, i);
                val2 = get_pixel(b, k, j, i);
                set_pixel(new_im, k, j, i, val1 + val2);
            }
        }
    }
    return new_im;
}

image sub_image(image a, image b)
{
    float val1, val2;
    if((a.h != b.h)&&(a.w != b.w)){
        return make_image(1,1,1);
    }
    image new_im;
    new_im = make_image(a.w, a.h, a.c);
    for(int i = 0; i<a.c; i++){
        for(int j = 0; j<a.h; j++){
            for(int k = 0; k<a.w; k++){
                val1 = get_pixel(a, k, j, i);
                val2 = get_pixel(b, k, j, i);
                set_pixel(new_im, k, j, i, val1 - val2);
            }
        }
    }
    return new_im;
}

image make_gx_filter()
{
    image im = make_image(3,3,1);
    im.data[0] = -1;
    im.data[1] = 0;
    im.data[2] = 1;
    im.data[3] = -2;
    im.data[4] = 0;
    im.data[5] = 2;
    im.data[6] = -1;
    im.data[7] = 0;
    im.data[8] = 1;
    return im;
}

image make_gy_filter()
{
    image im = make_image(3,3,1);
    im.data[0] = -1;
    im.data[1] = -2;
    im.data[2] = -1;
    im.data[3] = 0;
    im.data[4] = 0;
    im.data[5] = 0;
    im.data[6] = 1;
    im.data[7] = 2;
    im.data[8] = 1;
    return im;
}

void feature_ext(image im, float feature_val[]){
    float min_val = 1, max_val = 0;

    for(int i = 0; i < im.w; i++){
        for(int j = 0; j < im.h; j++){
            for(int k = 0; k < im.c; k++){
                if(get_pixel(im, i, j, k) < min_val) min_val = get_pixel(im, i, j, k);
                if(get_pixel(im, i, j, k) > max_val) max_val = get_pixel(im, i, j, k);
            }
        }
    }

    feature_val[0] = min_val;
    feature_val[1] = max_val - min_val;
}

void feature_normalize(image im)
{
    float feature_val[2];

    

    printf("%d, %d, %d\n", im.w, im.h, im.c);

    feature_ext(im, feature_val);
    float min_val = feature_val[0];
    float range = feature_val[1];
    float val;

    //

    for(int i = 0; i < im.w; i++){
        for(int j = 0; j < im.h; j++){
            for(int k = 0; k < im.c; k++){
                val = get_pixel(im, i, j, k);
                val = (val - min_val)/range;
                set_pixel(im, i, j, k, val);
            }
        }
    }

    
}

image *sobel_image(image im)
{
    image flat_im;
    image *res_img;
    image im1, im2, filter1, filter2;
    image Gx, Gy;

    float Gx_val, Gy_val;

    res_img = calloc(2, sizeof(image));

    im1 = make_image(im.w,im.h,im.c);
    im2 = make_image(im.w,im.h,im.c);

    Gx = copy_image(im);
    Gy = copy_image(im);

    filter1 = make_gx_filter();
    filter2 = make_gy_filter();

    Gx = convolve_image(Gx, filter1, 1);
    Gy = convolve_image(Gy, filter2, 1);
    
    for(int i = 0; i < im.w; i++){
        for(int j = 0; j < im.h; j++){
            for(int k = 0; k < im.c; k++){
                Gx_val = get_pixel(Gx, i, j, k);
                Gy_val = get_pixel(Gy, i, j, k);
                set_pixel(im1, i, j, k, sqrt(pow(Gx_val,2) + pow(Gy_val,2)));
                set_pixel(im2, i, j, k, atan2(Gx_val,Gy_val));
            }
        }
    }

    res_img[0] = im1;
    res_img[1] = im2;

    return res_img;
}

image colorize_sobel(image im)
{
    // TODO
    return make_image(1,1,1);
}
