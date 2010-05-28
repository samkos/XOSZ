#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

main(argc, argv)
int argc;
char *argv[];
{
    FILE *in;
    char flag[260],input[260],fic_controle[260];
    unsigned char c;

    for (c=0;c<255;c++) flag[c]=0;

    strcpy(fic_controle,argv[1]);
    
/*      printf("%d %s %s",argc,argv[0],argv[1]); */

   if (!(in=fopen(fic_controle,"r"))) 
     { printf ("\n pb ouverture fichier input "); for (c=1;c<255;) 
       exit(1);
     }
    fscanf(in,"%s",&input);
    fclose(in);
    for (c=0;c<strlen(input);c++) flag[input[c]]=1;

    for (flag['Q']=0;!flag['Q'];)
      { printf("\n --> ");
	for (c=0;c<255;c++) if (flag[c]&&c!='1') {printf("%c",c);};
	printf(" < \n ? -> ");
	scanf("%s",&input);
	for (c=0;c<strlen(input);c++) flag[input[c]]=1-flag[input[c]];
	flag['1']=1;
	if (!(in=fopen(fic_controle,"w"))) 
	  { printf ("\n pb ouverture fichier input ");
	    exit(1);
	  }
	for (c=0;c<255;c++) if (flag[c]) {fprintf(in,"%c",c);};
	fprintf(in,"\n\n");
	fclose(in);
      }

    return 0;
}



