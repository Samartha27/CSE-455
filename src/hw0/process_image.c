#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    
	if (x < 0) {x = 0;}
    	else if (x >= im.w){ x = im.w - 1;}
    	if (y < 0){ y = 0;}
    	else if (y >= im.h) {y = im.h - 1;}
	int index= c*im.h*im.w +y*im.w +x;
	float pixel_intensity =im.data[index];
	return pixel_intensity;
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
	if (x >= 0 && x < im.w && y >= 0 && y < im.h && c >= 0 && c < im.c){
	//printf("I am Here \n");
	int index = c*im.h*im.w + y*im.w + x;
	im.data[index] = v;
    }
}

image copy_image(image im)
{
	image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
	memcpy(copy.data, im.data,im.w*im.h*im.c* sizeof(float));		
	return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
     
    for (int y = 0; y < im.h; y++){
		for (int x = 0;x<im.w;x++){
			float r = get_pixel(im,x,y,0);
			float g = get_pixel(im,x,y,1);
			float b = get_pixel(im,x,y,2);
			float v = 0.299*r + 0.587*g + 0.114*b;
			set_pixel(gray,x,y,0,v);

		 }}
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int y = 0; y < im.h; y++){ 
	for (int x = 0;x<im.w;x++){
	    set_pixel(im,x,y,c,get_pixel(im,x,y,c)+v);}} 
}

void clamp_image(image im)
{
    // TODO Fill this in
    for (int c=0;c<im.c;c++){
	    for (int y =0;y<im.h;y++){
		    for (int x=0;x<im.w;x++){
			    float v = get_pixel(im,x,y,c);
			    if (v <0.0) {v = 0.0;}
			    else if (v >=1.0) {v =1.0;}
			    set_pixel(im,x,y,c,v);}}}
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    float H, S, V,red, green, blue,C, m = 0.0;

    for(int y = 0; y<im.h; y++){
        for(int x = 0; x<im.w; x++){
            red = get_pixel(im, x, y, 0);
	    	green = get_pixel(im, x, y, 1);
	    	blue = get_pixel(im, x, y, 2);
			H = 0;S = 0;V = 0;

	    	V = three_way_max(red, green, blue);
	    	m = three_way_min(red, green, blue);
	    	C = V - m;
	    
			if(V != 0){S = C/V;}
				if(C > 0){
					if(V == red){H = (green - blue)/C;}
					if(V == green){H = (blue - red)/C + 2;}
					if(V == blue){H = (red - green)/C + 4;}
					if(H < 0){H = H/6 + 1;}
					else{H = H/6;}}
			set_pixel(im, x, y, 0, H);
			set_pixel(im, x, y, 1, S);
			set_pixel(im, x, y, 2, V);
			
			//printf("H %f S %f V %f \n",H,S,V);
			//printf("R %f G %f B %f \n",red,green,blue);
	}
    }
}
void hsv_to_rgb(image im)
{
    // TODO Fill this in
    for(int y =0;y<im.h;y++){        

	    for (int x=0;x<im.w;x++){
		   float H =(get_pixel(im,x,y,0)*360.0);
	    	   float S = get_pixel(im,x,y,1);
		   float V = get_pixel(im,x,y,2);
		   float C = V*S;
		   H =( H/60.0);
		   float X= C*(1- fabs(fmod(H,2) -1));
		   float R = 0.0;
		   float G = 0.0;
		   float B = 0.0;
		   if (H>=0 &&H<1) {R=C;G=X;}
		   else if (H>=1 && H<2) {R=X; G=C;}
		   else if (H>=2 && H<3) {G=C;B=X;}
		   else if (H>=3 && H<4) {G=X;B=C;}
		   else if (H>=4 && H<5) {R=X;B=C;}
		   else {R=C;B=X;}
		   float m = V-C;
		   R=R+m;G=G+m;B=B+m;
		   set_pixel(im,x,y,0,R);
		   set_pixel(im,x,y,1,G);
		   set_pixel(im,x,y,2,B);
   		   //printf("H %f S %f V %f \n",H,S,V);
  		   //printf("R %f G %f B %f \n", R,G,B);

}}}

void scale_image(image im,int c,float v)
{
for (int y = 0; y < im.h; y++){                                                              for (int x = 0;x<im.w;x++){                                                                  set_pixel(im,x,y,c,get_pixel(im,x,y,c)*v);}}



}
