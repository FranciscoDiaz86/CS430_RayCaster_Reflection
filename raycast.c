#include <stdlib.h>
#include <math.h>
#include "csv_parser.h"
#include "3dmath.h"

typedef struct Pixel{
  unsigned char r, g, b;
}Pixel;

typedef struct PPMImage{
  int width, height;
  Pixel *pix_info;
}PPM_Image;

double* raycast (V3 Rd, V3 Ro, int obj_count, int light_count, Object *array, Object*light_array, double* background);
Pixel* render (int width, int height, Object camera, Object *obj_array, int obj_count, Object *light_array, int light_count);
double sphere_intersection(V3 Rd, V3 Ro, Object obj);
double plane_intersection(V3 Rd, V3 Ro, Object obj);
double* idiff(Object obj, Object light, double* intersection);
double* ispec(Object obj, Object light, double* intersection);
double clamp(double v);

int main (int argc, char *argv[]) {
  //checking for the right amout of arguments
  if (argc != 5) {
    fprintf(stderr, "Usage: The program needs a width height cvs file output ppm file.\n");
    return(1);
  }
  
  int obj_count = 0;
  int light_count = 0;
  Object light_array[128];
  Object obj_array[128];
  int width = atoi(argv[1]);
  int height = atoi(argv[2]);
  Pixel *pix_image;
  FILE *csv_file;
  FILE *output_file;
  
  csv_file = fopen(argv[3], "r");
  output_file = fopen(argv[4], "w");

  //Checking to see if the files opened without a problem
  if (csv_file == NULL) {
    fprintf(stderr, "Error: Problem opening the CSV file!\n");
    return(1);
  }

  if (output_file == NULL) {
    fprintf(stderr, "Error: Problem opening the output PPM file!\n");
    return(1);
  }

  fseek(csv_file, 0, SEEK_END);
  unsigned long len = (unsigned long)ftell(csv_file);
  if (len <= 0) {
    fprintf(stderr, "ERROR: The CSV file is empty!\n");
    return(1);
  }
  rewind(csv_file);
  
  while (1){
    if (feof(csv_file)){
	break;
    }
    
    Object obj;
    Object *obj_ptr;
    obj_ptr = &obj;
    parse_type(csv_file, obj_ptr);
    parse_field(csv_file, obj_ptr);
    obj_array[obj_count] = obj;
    obj_array[obj_count].pos = obj_count;
    obj_count++; 
  }

  for (int i = 0; i < obj_count; i++) {
    if (obj_array[i].kind == 4) {
      light_array[light_count] = obj_array[i];
      light_count++;
    }
  }
  
  if (light_count == 0) {
    fprintf(stderr, "ERROR: The file had no light objects in it!\n");
    return(1);
  }

  if (obj_count == 0) {
    fprintf(stderr, "ERROR: The file had no objects in it!\n");
    return(1);
  }
  
  Object camera_ptr;
  for (int i = 0; i < obj_count; i++) {
    if (obj_array[i].kind == 1) {
      camera_ptr = obj_array[i];
      break;
    }
    else {
      fprintf(stderr, "ERROR: The file had no camera object! Please add one.\n");
      return(1);
    }
  }


  pix_image = render(width, height, camera_ptr, obj_array, obj_count, light_array, light_count);
  
  //Writing P3 to file
  fprintf(output_file, "%s\n", "P3");

  //Writing width and height to file
  fprintf(output_file, "%s "  "%s\n", argv[1], argv[2]);

  //Writing max color value to file
  fprintf(output_file, "%s\n", "255");

  //Writing image data to file
  for (int i = 0; i < width * height; i++){
    fprintf(output_file, "%d\n", pix_image[i].r);
    fprintf(output_file, "%d\n", pix_image[i].g);
    fprintf(output_file, "%d\n", pix_image[i].b);
  }
  fclose(output_file);
  fclose(csv_file);
  
  return(0);
}
int raycast_primitive(V3 Rd, V3 Ro, Object closest_object){
  double t;
  if (closest_object.kind == 2){
    //sphere intersection
    t = sphere_intersection(Rd, Ro, closest_object);
  }

  if (closest_object.kind == 3){
    //plane intersection
    t = plane_intersection(Rd, Ro, closest_object);
  }

  if (t == INFINITY){
    return(0);
  }else{
    return(1);
  }
}

double* raycast (V3 Rd, V3 Ro, int obj_count, int light_count, Object *array, Object *light_array, double* background){
  Object closest_obj;
  double* background_color = background; 
  double* color_array = malloc(sizeof(double) * 3);
  double closest_t = INFINITY;
  double t;
  Object obj;
  
  for (int i = 0; i < obj_count; i++){
    obj = array[i];
    if (obj.kind == 2){
      //sphere intersection
      t = sphere_intersection(Rd, Ro, obj);
      if (t < closest_t){
	closest_t = t;
	closest_obj = array[i];
      }
    }
    
    if (obj.kind == 3){
      //plane intersection
      t = plane_intersection(Rd, Ro, obj);
      if (t < closest_t){
	closest_t = t;
	closest_obj = array[i];
      }
    }
  }

  //Modify here
  double* Rdt = malloc(sizeof(double) * 3);
  double* intersection_point = malloc(sizeof(double) * 3);
  double* Ro2 = malloc(sizeof(double) * 3);
  double* Rd2 = malloc(sizeof(double) * 3);
  double* Rds = malloc(sizeof(double) * 3);
  double* Ros = malloc(sizeof(double) * 3);
  double* Vo = malloc(sizeof(double) * 3);
  double* diff = malloc(sizeof(double) * 3);
  double* spec = malloc(sizeof(double) * 3);
  double ambient = 0.01;
  double emittence = 0.01;
  double* N = malloc(sizeof(double) * 3);
  double* E = malloc(sizeof(double) * 3);
  double f_rad, f_ang, x_inter, y_inter, z_inter;
  double t_light = INFINITY;
  int shadow = 0;
  Object closest_shadow_object;
  
  v3_assign(color_array, 0, 0, 0);

  if (closest_t != INFINITY){
    //Getting intersection point
    //v3_scale(Rdt, Rd, closest_t);
    //v3_add(intersection_point, Ro, Rdt);
    x_inter = Ro[0] + Rd[0] * closest_t;
    y_inter = Ro[1] + Rd[1] * closest_t;
    z_inter = Ro[2] + Rd[2] * closest_t;
    v3_assign(intersection_point, x_inter, y_inter, z_inter);
    
    //Calculating normal for sphere
    if (closest_obj.kind == 2){
      double x, y, z;
      x = (intersection_point[0] - obj.position[0])/obj.radius;
      y = (intersection_point[1] - obj.position[1])/obj.radius;
      z = (intersection_point[2] - obj.position[2])/obj.radius;
      v3_assign(N, x, y, z);
    }

    //Getting Normal for plane
    if (closest_obj.kind == 3){
      N = closest_obj.normal;
    }

    //Calculating epsilone value
    v3_assign(E, N[0]/2, N[1]/2, N[2]/2);
    
    for (int i = 0; i < light_count; i++){
      Object l = light_array[i];
      double dl = v3_distance(l.position, intersection_point);
      closest_shadow_object.is_null = 1;
      
      //Checking for shadows
      v3_sub(Rds, l.position, intersection_point);
      double smag = sqrt(pow(Rds[0], 2) + pow(Rds[1], 2) + pow(Rds[2], 2));
      v3_assign(Rds, Rds[0]/smag, Rds[1]/smag, Rds[2]/smag);
      v3_add(Ros, intersection_point, E);
      
      for(int j = 0; j < obj_count; j++){
	obj = array[j];

	if (array[j].pos == closest_obj.pos)
	  continue;
	
	if (obj.kind == 2){
	  //sphere intersection
	  t = sphere_intersection(Rds, Ros, obj);
	  if (t < t_light){
	    t_light = t;
	    closest_shadow_object = array[j];
	  }
	}
	
	if (obj.kind == 3){
	  //plane intersection
	  t = plane_intersection(Rds, Ros, obj);
	  if (t < t_light){
	    t_light = t;
	    closest_shadow_object = array[j];
	  }
	}
	if (t_light > dl)
	  continue;
      }
      
      if (closest_shadow_object.is_null == 1){
	//Adding epsilone value to Ro2
	v3_add(Ro2, intersection_point, E);
	
	//Normalizing Rd2
	v3_sub(Rd2, intersection_point, l.position);
	double mag = sqrt(pow(Rd2[0], 2) + pow(Rd2[1], 2) + pow(Rd2[2], 2));
	v3_assign(Rd2, Rd2[0]/mag, Rd2[1]/mag, Rd2[2]/mag);
	
	//checking to see if light hits object
	int hit = raycast_primitive(Rd2, Ro2, closest_obj);
	if (!hit) continue;
	f_rad = 1.0;
	f_ang = 1.0;
	
	v3_sub(Vo, intersection_point, l.position);
	// Radial attenuation
	if (dl != INFINITY) {
	  f_rad = 1 / (l.radial_a2 * dl * dl + l.radial_a1 * dl + l.radial_a0);
	}
	
	// Angular attenuation
	if (l.theta = 0){
	  f_ang = 1.0;
	}else if (v3_dot(l.direction, Vo) > cos(l.theta * (M_PI/180))){
	  f_ang = 0.0;
	}else{
	  f_ang = pow(v3_dot(Vo, l.direction), l.angular_a0);
	}
      
	//Adding the light values
	diff = idiff(closest_obj, l, intersection_point);
	spec = ispec(closest_obj, l, intersection_point);
	color_array[0] = f_rad * f_ang * (diff[0] + spec[0]);
	color_array[1] = f_rad * f_ang * (diff[1] + spec[1]);
	color_array[2] = f_rad * f_ang * (diff[2] + spec[2]);
      }
    }
    color_array[0] = color_array[0] + ambient + emittence;
    color_array[1] = color_array[1] + ambient + emittence;
    color_array[2] = color_array[2] + ambient + emittence;
    return(color_array);
  }else{
    return(color_array);
  }
}

double* idiff(Object obj, Object light, double* intersection){
  // Allocating space for light vector
  double* L = malloc(sizeof(double) * 3);
  double* diff_array = malloc(sizeof(double) * 3);
  double dot, r, g, b;
  v3_assign(diff_array, 0, 0, 0);
  // Calculating light vector
  v3_sub(L, light.position, intersection);
  
  if (obj.kind == 2){
    // Creating normal vector
    double* N = malloc(sizeof(double) * 3);
    double x, y, z;
    x = (intersection[0] - obj.position[0])/obj.radius;
    y = (intersection[1] - obj.position[1])/obj.radius;
    z = (intersection[2] - obj.position[2])/obj.radius;
    v3_assign(N, x, y, z);

    dot = v3_dot(N, L);
    if (dot > 0) {
      r = obj.diffuse[0] * light.color[0] * dot;
      g = obj.diffuse[1] * light.color[0] * dot;
      b = obj.diffuse[2] * light.color[2] * dot;

      v3_assign(diff_array, r, g, b);
      return(diff_array);
    }else{
      v3_assign(diff_array, 0, 0, 0);
      return(diff_array);
    }
    return(diff_array);
  }
 
  if (obj.kind == 3) {
    dot = v3_dot(obj.normal, L);
    if (dot > 0) {
      r = obj.diffuse[0] * light.color[0] * dot;
      g = obj.diffuse[1] * light.color[0] * dot;
      b = obj.diffuse[2] * light.color[2] * dot;

      v3_assign(diff_array, r, g, b);
      return(diff_array);
    }else{
      v3_assign(diff_array, 0, 0, 0);
      return(diff_array);
    }
  }
}

double* ispec(Object obj, Object light, double* intersection){
  double* V = malloc(sizeof(double) * 3);
  double* R = malloc(sizeof(double) * 3);
  double* spec_array = malloc(sizeof(double) * 3);
  double vn, r, g, b, dot;
  v3_sub(V, light.position, intersection);

  if (obj.kind == 2){
    //Creating the normal
    double* N = malloc(sizeof(double) * 3);
    double x, y, z;
    x = (intersection[0] - obj.position[0])/obj.radius;
    y = (intersection[1] - obj.position[1])/obj.radius;
    z = (intersection[2] - obj.position[2])/obj.radius;
    v3_assign(N, x, y, z);

    //Normalizing v
    double mag = sqrt(pow(V[0], 2) + pow(V[1], 2) + pow(V[2], 2));
    v3_assign(V, V[0]/mag, V[1]/mag, V[2]/mag);
    
    //Creating R
    vn = v3_dot(V, N);
    vn = vn * 2;
    v3_scale(N, N, vn);
    v3_sub(R, N, V);

    dot = v3_dot(V, R);
    if (dot > 0){
      dot = pow(dot, 20);
      r = obj.specular[0] * light.color[0] * dot;
      g = obj.specular[1] * light.color[1] * dot;
      b = obj.specular[2] * light.color[2] * dot;

      v3_assign(spec_array, r, g, b);
      return(spec_array);
    }else{
      v3_assign(spec_array, 0, 0, 0);
      return(spec_array);
    }
  }

  if (obj.kind == 3){
    v3_assign(spec_array, 0, 0, 0);
    return(spec_array);   
  }
}
    
    
    
Pixel* render (int width, int height, Object camera, Object *obj_array, int obj_count, Object *light_array, int light_count){
  Pixel *pix_info;
  double *color;
  double background[3] = {0.97647, 0.513725, 0.12549};
  double dx = (camera.width/width);
  double dy = (camera.height/height);
  int current_pix = 0;
  
  pix_info = (Pixel *)malloc(sizeof(Pixel) * width * height);
  double* Pij = malloc(sizeof(double) * 3);
  double* Rd = malloc(sizeof(double) * 3);
  double* Ro = malloc(sizeof(double) * 3);
  v3_assign(Ro, 0 , 0, 0);
  for (int i = 0; i < height; i++){
    double Py = 0 - (camera.height/2) + dy * (i + 0.5);
    for (int j = 0; j < width; j++){
      double Px = 0 - (camera.width/2) + dx * (j + 0.5);
     
      v3_assign(Pij,
		Px,
		-Py,
		-1);
      
      v3_sub(Rd, Pij, Ro);
      
      double mag = sqrt(pow(Rd[0], 2) + pow(Rd[1], 2) + pow(Rd[2], 2));

     
      v3_assign(Rd,
		Rd[0]/mag,
		Rd[1]/mag,
		Rd[2]/mag);

      color = raycast(Rd, Ro, obj_count, light_count, obj_array, light_array, background);
      pix_info[current_pix].r = (clamp(color[0]) * 255);
      pix_info[current_pix].g = (clamp(color[1]) * 255);
      pix_info[current_pix].b = (clamp(color[2]) * 255);
      current_pix++;
    }
  }
  return(pix_info);
}

double clamp(double v){
  if (v > 1.0)
    return(1.0);
  else if (v < 0.0)
    return(0.0);
  else
    return(v);
}

double sphere_intersection (V3 Rd, V3 Ro, Object obj){
  //Rd values for x y z
  double Rdx = Rd[0];
  double Rdy = Rd[1];
  double Rdz = Rd[2];

  //R0 values
  double Rox = Ro[0];
  double Roy = Ro[1];
  double Roz = Ro[2];
  
  //Sphere values for x y z and r
  double Cx = obj.position[0];
  double Cy = obj.position[1];
  double Cz = obj.position[2];
  double r = obj.radius;
  //Calculating A B C for quadratic equation
  double A = pow(Rdx, 2) + pow(Rdy, 2) + pow(Rdz, 2);
  double B = 2 * (Rdx * (Rox - Cx) + Rdy * (Roy - Cy) + Rdz * (Roz - Cz));
  double C = ((Rox - Cx) * (Rox - Cx)) + ((Roy - Cy) * (Roy - Cy)) + ((Roz - Cz) * (Roz - Cz)) - (r * r);

  //Initializing t0 t1
  double t0;
  double t1;

  //Calculating discriminant
  double dis;
  dis = pow(B,2) - (4 * C); 
  if (dis < 0){
    return(INFINITY);
  }
  
  //Calculating for t0 and t1
  double sqrt_dis = sqrt(dis);
  t0 = (-B - sqrt_dis)/2;
  if (t0 > 0){
    return(t0);
  }else{
    t1 = (-B + sqrt(dis))/2;
    return(t1);
  }
}

double plane_intersection (V3 Rd, V3 Ro, Object obj){
  double Vd, Vo, D, t;
  D = v3_distance(Ro, obj.position);
  Vd = v3_dot(obj.normal, Rd);
  if (Vd < 0){
    Vo = -(v3_dot(obj.normal, Ro) + D);
    t = Vo/Vd;
    if (t < 0)
      return INFINITY;
    else
      return t;
  }
  return INFINITY;
}
  
