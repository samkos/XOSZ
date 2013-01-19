#include <string.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <unistd.h>

#define PANIC(a) do { \
                perror(a); \
                if (temp_name) unlink(temp_name);\
                exit(1);\
        } while(0)

#define create 1
#define close  5

#if SGI
#if petit
#define geom "gnuplot -geometry 172x109+1108+%d"
#else
#define geom "gnuplot -geometry 344x218+0936+%d"
#define vim  "xterm -geometry 10x5+730+0 -title '%s : zef_mon' -e ./zplot_control %s"
#define add_geom 289
#endif
#define minnew 0
#else
#if petit
#define vim  "xterm -geometry 10x5+534+0 -title '%s : zef_mon' -e ./zplot_control %s"
#define geom "gnuplot -geometry 172x109+608+%d"
#define add_geom 138
#else
#define vim  "xterm -geometry 10x5+421+00 -title '%s : zef_mon' -e ./zplot_control %s"
#define geom "gnuplot -geometry 300x156+496+%d"
#define add_geom 185
#endif
#define minnew 15
#endif 

#define one_d   1
#define contour 2
#define surface 10

int plot_ (what,length_what,splot,valmin,valmax) 
    int *length_what, *splot;
    char *what;
    double *valmin,*valmax;
{
    static FILE *command[260];
    static FILE *controle,*vi;
    double a,b;
    int i,j,k;
    char s[100],c,texte[100],s2[100];
    static int new=0,deja[260];


    if (!new&&*splot==100) return 0;

    if (!new)
      { for (i=0;i<259;i++) deja[i]=0;
        new=1;
	for (i=0;i<*length_what;i++)
	  { s[i]=what[i];
	    if (s[i]=='/') {j=i;}
	    if (s[i]=='.') {c= *(s+i-1); k=i;}
	  }
	s[i]=0;

	s[k-1]=0;
	strcpy(s2,s+j+1);
	s[k-1]='Z';
	sprintf(texte,vim,s2,s);
	vi = popen(texte,"w"); 

	s[k-1]=c;
	
     }
      

    if (*splot==100) 
      { c=what[0];
	if (deja[c])  { fclose(command[c]); new-=add_geom;  }
	deja[c]=0;
	return 0;
      }
    


    for (i=0;i<*length_what;i++)
      { s[i]=what[i];
        if (s[i]=='/') {j=i;}
        if (s[i]=='.') {c= *(s+i-1);}
      }
    s[i]=0;

    if (!deja[c]) {
      sprintf(texte,geom,new-minnew);
      new += add_geom; deja[c]=1;
      command[c] = popen(texte,"w"); 
      switch (*splot) {
      case (contour) : fprintf(command[c],"set view 0,0,1.7,1.7\n set nosurface\n set cntrparam levels 28 \
                           \n set contour\n set nokey\n set title \"%s\" \
                           \n splot '%s' w line\n ",s+j+1,s); 
                       break;
      case (surface) : fprintf(command[c],"set contour\n set nokey\n set title \"%s\" \
                           \n splot '%s' w line\n ",s+j+1,s); 
                       break;
      case (one_d) :   fprintf(command[c],"set nokey\n set title \"%s\" \
                           \n set view 90,0,1.7,1. \n splot '%s' w line\n ",s+j+1,s); 
                       break;
      
      }
    }
    else
      { fprintf(command[c],"replot \n "); }
    
    fflush(command[c]);

    return 0;
}


