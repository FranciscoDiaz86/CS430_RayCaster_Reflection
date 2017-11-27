#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct{
  int kind, is_null, pos;
  double position[3], color[3], normal[3], diffuse[3], specular[3], direction[3];
  double width, height, radius, theta, angular_a0, radial_a0, radial_a1, radial_a2, reflect, refrac, ior;
}Object;

int parse_width(FILE *fh, Object *obj);
int parse_height(FILE *fh, Object *obj);
int parse_color(FILE *fh, Object *obj);
int parse_position(FILE *fh, Object *obj);
int parse_normal(FILE *fh, Object *obj);
int parse_radius(FILE *fh, Object *obj);
int parse_diffuse(FILE *fh, Object *obj);
int parse_specular(FILE *fh, Object *obj);
int parse_theta(FILE *fh, Object *obj);
int parse_radial_a0(FILE *fh, Object *obj);
int parse_radial_a1(FILE *fh, Object *obj);
int parse_radial_a2(FILE *fh, Object *obj);
int parse_direction(FILE *fh, Object *obj);
int parse_angular_a0(FILE *fh, Object *obj);
int parse_reflectivity(FILE *fh, Object *obj);
int parse_refractivity(FILE *fh, Object *obj);
int parse_ior(FILE *fh, Object *obj);

int parse_type(FILE *fh, Object *obj){
  int chr;
  char buff[255];
  int buff_pos = 0;
  int file_ptr;
   
  while(1){
    chr = fgetc(fh);
      
    if (chr == ','){
      break;
    }

    buff[buff_pos] = chr;
    buff[buff_pos + 1] = '\0';
    buff_pos++;
  }
  
  if((strcmp(buff, "camera")) == 0){
    obj -> kind = 1;
  }

  if((strcmp(buff, "sphere")) == 0){
    obj -> kind = 2;
  }

  if((strcmp(buff, "plane")) == 0){
    obj -> kind = 3;
  }
  
  if((strcmp(buff, "light")) == 0){
      obj -> kind = 4;
  }
  
  return(0);
}

int parse_field(FILE *fh, Object *obj){
  int chr;
  char buff[255];
  int buff_pos = 0;

  chr = fgetc(fh);
  while(chr != EOF){
    
    if(chr == ' '){
      chr = fgetc(fh);
    }
    
    if(chr != ' ' && chr != ':'){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
    
    if(chr == ':'){
      if((strcmp(buff, "width")) == 0){
	int current_fptr = parse_width(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
	  
	  if((strcmp(buff, "reflectivity")) == 0){
	int current_fptr = parse_reflectivity(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
	  
	  if((strcmp(buff, "refractivity")) == 0){
	int current_fptr = parse_refractivity(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
	  
	  if((strcmp(buff, "ior")) == 0){
	int current_fptr = parse_ior(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "height")) == 0){
	int current_fptr = parse_height(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "color")) == 0){
	int current_fptr = parse_color(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "position")) == 0){
	int current_fptr = parse_position(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }

      if((strcmp(buff, "radius")) == 0){
	int current_fptr = parse_radius(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }

      if((strcmp(buff, "normal")) == 0){
	int current_fptr = parse_normal(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }

      if((strcmp(buff, "direction")) == 0){
	int current_fptr = parse_direction(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "diffuse_color")) == 0){
	int current_fptr = parse_diffuse(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "specular_color")) == 0){
	int current_fptr = parse_specular(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "theta")) == 0){
	int current_fptr = parse_theta(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "radial-a0")) == 0){
	int current_fptr = parse_radial_a0(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "radial-a1")) == 0){
	int current_fptr = parse_radial_a1(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
      
      if((strcmp(buff, "radial-a2")) == 0){
	int current_fptr = parse_radial_a2(fh, obj);
	if(current_fptr == '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }

      if((strcmp(buff, "angular-a0")) == 0){
	int current_fptr = parse_angular_a0(fh, obj);
	if(current_fptr = '\n'){
	  break;
	}
	chr = fgetc(fh);
	buff_pos = 0;
      }
    }
  }
  return(0);
}

int parse_width(FILE *fh, Object *obj){
  double obj_width;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_width = atof(buff);
  obj -> width = obj_width;
  return(chr);
}

int parse_reflectivity(FILE *fh, Object *obj){
  double obj_reflect;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_reflect = atof(buff);
  obj -> reflect = obj_reflect;
  return(chr);
}

int parse_refractivity(FILE *fh, Object *obj){
  double obj_refract;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_refract = atof(buff);
  obj -> refrac = obj_refract;
  return(chr);
}

int parse_ior(FILE *fh, Object *obj){
  double obj_ior;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_ior = atof(buff);
  obj -> ior = obj_ior;
  return(chr);
}

int parse_height(FILE *fh, Object *obj){
  double obj_height;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_height = atof(buff);
  obj -> height = obj_height;
  return(chr);
}

int parse_color(FILE *fh, Object *obj){
  double obj_color;
  char buff[255];
  int buff_pos = 0;
  int clr_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }

    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr == '['){
      chr = fgetc(fh);
      while(1){
	if(chr == ']'){
	  obj_color = atof(buff);
	  obj -> color[clr_pos] = obj_color;
	  chr = fgetc(fh);
	  break;
	}
	if(chr == ' '){
	  chr = fgetc(fh);
	}

	if(chr == ','){
	  obj_color = atof(buff);
	  obj -> color[clr_pos] = obj_color;
	  clr_pos++;
	  buff_pos = 0;
	  chr = fgetc(fh);
	}

	if(chr != ' ' && chr != ','){
	  buff[buff_pos] = chr;
	  buff[buff_pos + 1] = '\0';
	  buff_pos++;
	  chr = fgetc(fh);
	}
      }
      break;
    }
  }
  chr = fgetc(fh);
  return(chr);
}

int parse_position(FILE *fh, Object *obj){
  double obj_position;
  char buff[255];
  int buff_pos = 0;
  int pos_ptr = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }

    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr == '['){
      chr = fgetc(fh);
      while(1){
	if(chr == ']'){
	  obj_position = atof(buff);
	  obj -> position[pos_ptr] = obj_position;
	  chr = fgetc(fh);
	  break;
	}
	if(chr == ' '){
	  chr = fgetc(fh);
	}

	if(chr == ','){
	  obj_position = atof(buff);
	  obj -> position[pos_ptr] = obj_position;
	  pos_ptr++;
	  buff_pos = 0;
	  chr = fgetc(fh);
	}

	if(chr != ' ' && chr != ','){
	  buff[buff_pos] = chr;
	  buff[buff_pos + 1] = '\0';
	  buff_pos++;
	  chr = fgetc(fh);
	}
      }
      break;
    }
  }
  chr = fgetc(fh);
  return(chr);
}

int parse_normal(FILE *fh, Object *obj){
  double obj_normal;
  char buff[255];
  int buff_pos = 0;
  int norm_ptr = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }

    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr == '['){
      chr = fgetc(fh);
      while(1){
	if(chr == ']'){
	  obj_normal = atof(buff);
	  obj -> normal[norm_ptr] = obj_normal;
	  chr = fgetc(fh);
	  break;
	}
	if(chr == ' '){
	  chr = fgetc(fh);
	}

	if(chr == ','){
	  obj_normal = atof(buff);
	  obj -> normal[norm_ptr] = obj_normal;
	  norm_ptr++;
	  buff_pos = 0;
	  chr = fgetc(fh);
	}

	if(chr != ' ' && chr != ','){
	  buff[buff_pos] = chr;
	  buff[buff_pos + 1] = '\0';
	  buff_pos++;
	  chr = fgetc(fh);
	}
      }
      break;
    }
  }
  chr = fgetc(fh);
  return(chr);
}

int parse_radius(FILE *fh, Object *obj){
  double obj_radius;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_radius = atof(buff);
  obj -> radius = obj_radius;
  return(chr);
}

int parse_diffuse(FILE *fh, Object *obj){
  double obj_diffuse;
  char buff[255];
  int buff_pos = 0;
  int diff_ptr = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }

    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr == '['){
      chr = fgetc(fh);
      while(1){
	if(chr == ']'){
	  obj_diffuse = atof(buff);
	  obj -> diffuse[diff_ptr] = obj_diffuse;
	  chr = fgetc(fh);
	  break;
	}
	if(chr == ' '){
	  chr = fgetc(fh);
	}

	if(chr == ','){
	  obj_diffuse = atof(buff);
	  obj -> diffuse[diff_ptr] = obj_diffuse;
	  diff_ptr++;
	  buff_pos = 0;
	  chr = fgetc(fh);
	}

	if(chr != ' ' && chr != ','){
	  buff[buff_pos] = chr;
	  buff[buff_pos + 1] = '\0';
	  buff_pos++;
	  chr = fgetc(fh);
	}
      }
      break;
    }
  }
  chr = fgetc(fh);
  return(chr);
}

int parse_specular(FILE *fh, Object *obj){
  double obj_specular;
  char buff[255];
  int buff_pos = 0;
  int spec_ptr = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }

    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr == '['){
      chr = fgetc(fh);
      while(1){
	if(chr == ']'){
	  obj_specular = atof(buff);
	  obj -> specular[spec_ptr] = obj_specular;
	  chr = fgetc(fh);
	  break;
	}
	if(chr == ' '){
	  chr = fgetc(fh);
	}

	if(chr == ','){
	  obj_specular = atof(buff);
	  obj -> specular[spec_ptr] = obj_specular;
	  spec_ptr++;
	  buff_pos = 0;
	  chr = fgetc(fh);
	}

	if(chr != ' ' && chr != ','){
	  buff[buff_pos] = chr;
	  buff[buff_pos + 1] = '\0';
	  buff_pos++;
	  chr = fgetc(fh);
	}
      }
      break;
    }
  }
  chr = fgetc(fh);
  return(chr);
}

int parse_theta(FILE *fh, Object *obj){
  double obj_theta;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_theta = atof(buff);
  obj -> theta = obj_theta;
  return(chr);
}

int parse_radial_a0(FILE *fh, Object *obj){
  double obj_a0;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_a0 = atof(buff);
  obj -> radial_a0 = obj_a0;
  return(chr);
}

int parse_radial_a1(FILE *fh, Object *obj){
  double obj_a1;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_a1 = atof(buff);
  obj -> radial_a1 = obj_a1;
  return(chr);
}

int parse_radial_a2(FILE *fh, Object *obj){
  double obj_a2;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_a2 = atof(buff);
  obj -> radial_a2 = obj_a2;
  return(chr);
}

int parse_angular_a0(FILE *fh, Object *obj){
  double obj_a0;
  char buff[255];
  int buff_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }
    
    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr != ' ' && chr != ','){
      buff[buff_pos] = chr;
      buff[buff_pos + 1] = '\0';
      buff_pos++;
      chr = fgetc(fh);
    }
  }
  obj_a0 = atof(buff);
  obj -> angular_a0 = obj_a0;
  return(chr);
}

int parse_direction(FILE *fh, Object *obj){
  double obj_direction;
  char buff[255];
  int buff_pos = 0;
  int clr_pos = 0;
  int chr;

  chr = fgetc(fh);
  while(chr != EOF){
    if(chr == ',' || chr == '\n'){
      break;
    }

    if(chr == ' '){
      chr = fgetc(fh);
    }

    if(chr == '['){
      chr = fgetc(fh);
      while(1){
	if(chr == ']'){
	  obj_direction = atof(buff);
	  obj -> direction[clr_pos] = obj_direction;
	  chr = fgetc(fh);
	  break;
	}
	if(chr == ' '){
	  chr = fgetc(fh);
	}

	if(chr == ','){
	  obj_direction = atof(buff);
	  obj -> direction[clr_pos] = obj_direction;
	  clr_pos++;
	  buff_pos = 0;
	  chr = fgetc(fh);
	}

	if(chr != ' ' && chr != ','){
	  buff[buff_pos] = chr;
	  buff[buff_pos + 1] = '\0';
	  buff_pos++;
	  chr = fgetc(fh);
	}
      }
      break;
    }
  }
  chr = fgetc(fh);
  return(chr);
}
