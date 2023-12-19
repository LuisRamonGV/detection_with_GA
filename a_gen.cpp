#include <cstdio>
#include <random>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <string>
#include <iostream>

//* ##########################################
//* Estructura del Algoritmo Genetico 
//* ##########################################
typedef struct{
    unsigned char   *Crom;  //* Cromosoma 
    int             *Vent;  //* Valores decodificados int
    double          Vfit;
    double          Vobj;

    double cx;  //* Centro X
    double cy;  //* Centro Y
    double r;   //* radio
} INDIVIDUO;

class GA{
private:
    INDIVIDUO          *Pob;       //* Poblacion int
    INDIVIDUO          *NewPob;    //* Poblacion int
    unsigned int        PobSize;   //* Tama;o de la poblacion
    unsigned int        NGens;     //* Numero de Genes
    const unsigned int *BitxGens;  //* Numero de Bits de cada gen
    unsigned int        CromSize;  //* Tama;o total de cromosoma
    double              SumaObj;
    double              SumaFit;   //* Suma total de 
    double              PromedioObj;
    unsigned int        IdBest;
    unsigned int        IdWorst;
    unsigned int        *Seleccionado;
    unsigned int        MaxGeneraciones;

    //* Procesamiento de la imagen
    unsigned int    NPix;
	int             AnchoIMG;
	int             AltoIMG;
	int             *Coorx;
	int             *Coory;

public:
    GA( const unsigned int NumGeneration, unsigned int SizePob, unsigned int NumDeGens, const unsigned int *NBitsXGen,
        int pix, int *c__x, int *c__y, int i_height, int i_width);
    ~GA();

    void ImprimeInd( unsigned int id );
    void imprimeResumen();
    void ImprimePoblacion();
    void DecodeEntero();
    void Obj_a_Fit();
	void Seleccion();
    void Cruza(double pC);
    void Muta(double pM);
    void Elitismo();
    void EvaluarPob();
    void NextGeneration();
    void process(double cruza, double muta);
    void ImprimeBest();

    double FuncionObjetivo(unsigned int i);

    void calcular_centro_radio(int k);
    void calcular_coeficientes(double &A, double &B, double &C, double &D, int Xi, int Yi, int Xj, int Yj, int Xk, int Yk, double &cx, double &cy);
    double calcular_distancia(double x1, double y1, double x2, double y2);
    void centro_radio();

};

//*                               i  j  k
int     centro_radio_values[3] = {0, 0, 0};
unsigned int NumPixNegros = 0, height, width;

//Para guardar el vector de coordenadas
int*    c_x; 
int*    c_y; 

GA::GA( const unsigned int NumGeneration, unsigned int SizePob, unsigned int NumDeGens, const unsigned int *NBitsXGen, 
        int pix, int* c__x, int* c__y, int i_height, int i_width){
    // Num. Generaciones / Tamano poblacion / Num. Genes / Num. bitsxGen / Num. Pixeles / coord. x / coord. y / alto / ancho
    unsigned int k, i;

    MaxGeneraciones = NumGeneration;
    PobSize  = SizePob;
    NGens    = NumDeGens;
    BitxGens = NBitsXGen;

    CromSize = 0;
    SumaObj = 0;
    PromedioObj = 0;
    IdBest = 0;
    IdWorst = 0;

    NPix        = pix;
    AltoIMG     = i_height;
    AnchoIMG    = i_width;
    Coorx      = c__x;
    Coory      = c__y;

    //* El tamaño total del cromosoma dependera del numero de bits
    for (k = 0; k<NGens; k++) CromSize+=BitxGens[k];

    //* Genera la memoria para una poblacion de N individuos
    Pob = new INDIVIDUO[PobSize];

    //* Generar la memoria para cada individuo
    for (k = 0; k<PobSize; k++) {
        Pob[k].Crom  = new unsigned char[CromSize];
        Pob[k].Vent  = new int[NGens];

        Pob[k].cx = 0;
        Pob[k].cy = 0;
        Pob[k].r  = 0;
    }
    
    //* Generar la memoria para la lista de Seleccionados
    Seleccionado = new unsigned int [PobSize];
    
    //* Generar la memoria para la Nueva Poblacion de N individuos
    NewPob = new INDIVIDUO[PobSize];

    //* Generar la memoria para cada individuo
    for (k = 0; k<PobSize; k++) {
        NewPob[k].Crom  = new unsigned char[CromSize];
        NewPob[k].Vent  = new int[NGens];
    }

    //! Inicializar el cromosoma
    for (k = 0; k<PobSize;  k++)    //* Para toda la poblacion
    for (i = 0; i<CromSize; i++)    //* Para todo el cromosoma del k-esimo individuo
        Pob[k].Crom[i] = rand()%2;
};

void GA::ImprimeInd( unsigned int id ){
    unsigned int k, acumulador, gen;
    gen = 0;
    acumulador = BitxGens[0];

    printf("\n%i: ", id);   
    //* Imprimir el cromosoma del individuo 
    for ( k = 0; k < CromSize; k++ ){
        if(k==acumulador){
            printf(":");
            gen++;
            acumulador += BitxGens[gen];
        }
        printf("%i", Pob[id].Crom[k]);
    }

    printf(" :\n");
    for ( k = 0; k < NGens; k++ ) printf(" \t%i", Pob[id].Vent[k]);
    
    printf(" | Vfit: %2.3f", Pob[id].Vfit);
    printf("\n--------------------------------------------------------------");
}

void GA::imprimeResumen(){
    printf("\nCentro: (%.2f,%.2f)", Pob[IdBest].cx, Pob[IdBest].cy);
	printf("\nRadio: %.2f", Pob[IdBest].r);
}

void GA::ImprimePoblacion(){
    for( unsigned int k = 0; k < PobSize; k++)
        ImprimeInd(k);
}

void GA::ImprimeBest(){
    ImprimeInd(IdBest);
};

void GA::DecodeEntero(){
    unsigned int k, i, j, acumulador, pos;

    for(k = 0; k<PobSize; k++){          //* Para toda la poblacion 
        j = 0;
        acumulador = 0;

        for(i = 0; i<NGens; i++){           //* Para cada Gen del k-esimo indice 
            Pob[k].Vent[i]=0;
            acumulador += BitxGens[i];
            pos = 0;

            for( ; j<acumulador; j++) { //* Para cada bit de cada gen
                Pob[k].Vent[i] += Pob[k].Crom[j]*pow(2,pos);
                pos++;
            }

            //Si el número entero es mayor al número de pixeles, asignarle un valor en el centro de la lista los pixeles 
            //* Revisar esto
            if(Pob[k].Vent[i] >= NPix)
                Pob[k].Vent[i] = NPix/2;
        }
    }
}
//________________________________________________________________________________________________________________________________________
double GA::FuncionObjetivo(unsigned int id) {
	    double fit=0;
		double xi, yi;
		int Ns = 40;

		for(double k=0; k<Ns;k++)
		{
			xi = (Pob[id].cx + (Pob[id].r*cos((k*M_PI)/Ns)));
			yi = (Pob[id].cy + (Pob[id].r*sin((k*M_PI)/Ns)));
			
            //* Redondeamos 
		    int Xint, Yint;
			if((xi-(int)xi)>=0.5)   {Xint = (int)xi + 1;}
			else                    {Xint = (int)xi;}

			if((yi-(int)yi)>=0.5)   {Yint = (int)yi + 1;}
			else                    {Yint = (int)yi;}\

			for(int j=0; j<NPix;j++)
				if((Xint==Coorx[j])&&(Yint==Coory[j])) fit++;
		}
		if((Pob[id].cx==0)&&(Pob[id].cy==0)) return 0;

		return fit/Ns;
	}

void GA::Obj_a_Fit(){
    unsigned int k;
    double rango, porcentaje;
    rango = Pob[IdBest].Vobj - Pob[IdWorst].Vobj;
    if(fabs(rango)>0.000001){
        SumaFit = 0;
        for (k = 0; k < PobSize; k++) {// Para toda la poblacion
            porcentaje = (Pob[k].Vobj - Pob[IdWorst].Vobj)/ rango;
            Pob[k].Vfit = 100*porcentaje;
            SumaFit += Pob[k].Vfit;
        }   
    } else {
        SumaFit = PobSize*100;
        for (k = 0; k < PobSize; k++) // Para toda la poblacion
            Pob[k].Vfit = 100;
    }
}

void GA::Seleccion(){
    double *Ruleta = new double [PobSize];
    double suma, pelota;
    int k, n;

    //* Crear la ruleta
    for(k = 0; k<PobSize; k++){ //* Para toda la poblacion 
        Ruleta[k] = Pob[k].Vfit/SumaFit;
        suma += Ruleta[k];
    }

    for (n = 0; n<PobSize; n++){
        //* Aventar la pelota
        pelota = (double)rand()/RAND_MAX;
        suma = 0;

        //* Calcular la seccion donde cae la pelota 
        for(k = 0; k<PobSize; k++){
            suma += Ruleta[k];
            if(pelota < suma){
                Seleccionado[n] = k;
                break;
            }
        }
    }
}

void GA::Cruza(double pC){
    int k, j, PdC, b;
    double r;
    for(k = 0; k<PobSize; k += 2){
        r = (double)rand()/RAND_MAX;
        if(r <= pC){ //* Se aplica la cruza
            PdC = rand()%CromSize;

            for(b = 0; b<=PdC; b++){
                NewPob[k].Crom[b]   = Pob[Seleccionado[k]].Crom[b];
                NewPob[k+1].Crom[b] = Pob[Seleccionado[k+1]].Crom[b];
            }
            for(b=PdC+1; b < CromSize; b++){
                NewPob[k+1].Crom[b] = Pob[Seleccionado[k]].Crom[b];
                NewPob[k].Crom[b]   = Pob[Seleccionado[k+1]].Crom[b];
            }

        } else {    //* No se aplica la cruza 
            for(j = 0; j<CromSize; j++){
                NewPob[k].Crom[j]   = Pob[Seleccionado[k]].Crom[j];
                NewPob[k+1].Crom[j] = Pob[Seleccionado[k+1]].Crom[j];
            }
        }
    }
}

void GA::Muta(double pM){
    int k, j, PdC, b;
    double r;
    for(k = 0; k<PobSize; k++)
    for(b = 0; b<CromSize; b++){
        r = (double)rand()/RAND_MAX;

        if(r <= pM){
            if(NewPob[k].Crom[b])   NewPob[k].Crom[b] = 0;
            else                    NewPob[k].Crom[b] = 1;
        }
    }
}

void GA::Elitismo(){
    unsigned int i;
    for(i = 0; i<CromSize; i++)
        NewPob[0].Crom[i]   = Pob[IdBest].Crom[i];
}

void GA::EvaluarPob(){
    double aux, BestObj, WorstObj;
    unsigned int k;
    // Evaluar el primer indiviuo o converit el 
    // valor objetivo en el mejor y el peor
    SumaObj = 0;
    IdWorst = 0;
    IdBest = 0;
    aux = FuncionObjetivo(0);
    Pob[0].Vobj = aux;
    BestObj = aux;
    WorstObj = aux;
    SumaObj += aux;  
    IdBest = 0;
    //SumaObj = 0;
    IdWorst = 0;
    for (k = 1; k < PobSize; k++){ // Para toda la poblacion
        aux = FuncionObjetivo(k);
        Pob[k].Vobj = aux;
        SumaObj += aux;
        if (aux > BestObj ){
            BestObj = aux;
            IdBest = k;
        }
        if (aux < WorstObj ){
            WorstObj = aux;
            IdWorst = k;
        }
    }
    PromedioObj = SumaObj / PobSize;
}

void GA::NextGeneration(){
    INDIVIDUO *Aux;
    Aux = Pob;

    Pob = NewPob;
    NewPob = Aux;
}

GA::~GA () {
    unsigned int k;

    //* Libera la memoria para cada individuo
    for(k = 0; k<PobSize; k++){
        delete[] Pob[k].Crom;
        delete[] Pob[k].Vent;
    }

    //* Libera la memoria de la poblacion 
    delete[] Pob;
};

void GA::process(double cruza, double muta){
    unsigned int t = 1;

    DecodeEntero();
    centro_radio();
    EvaluarPob();
    Obj_a_Fit();

    while(t<=MaxGeneraciones){
        Seleccion();
        Cruza(cruza);
        Muta(muta);
        Elitismo();

        NextGeneration();

        DecodeEntero();
        centro_radio();
        EvaluarPob();
        Obj_a_Fit();

        t++;
    }

    printf("\nGeneracion: %i", t-1);
    ImprimeBest();
    imprimeResumen();
}

//* ##########################################
//* PROCESAMIENTO DE IMAGENES
//* ##########################################
typedef unsigned char  		byte;    // Tipo de dato de 1 byte
typedef unsigned short int 	word;    // Tipo de dato de 2 bytes
typedef unsigned long  int 	dword;   // Tipo de dato de 4 bytes

typedef struct{ 
    byte    id[2];  // Identificador de fila BMP
    word    offset; // Offset al principio de la imagen
    word    ancho;  // Columnas de la imagen
    word    alto;   // Filas de la imagen
    byte    bpp;    // Bits de color por pixel
    int     size;   // Tamaño de la imagen
    byte    *head;  // Puntero al encabezado
    float   *imx;   // Puntero al inicio de la imagen
    } gcIMG;

gcIMG* 	    gcGetImgBmp(char *ruta);
void 	    gcPutImgBmp(char *ruta, gcIMG *img);
gcIMG*	    gcNewImg(int ancho,int alto);
void 	     gcFreeImg (gcIMG *img);

int contarBits(int numero) {
    int bits = 0;
    while (numero > 0) {
        bits++;
        numero = numero / 2;
    }
    return bits;
}  

void read(){
    //* Lee imagen
    //Declarar un Puntero a imagen
	gcIMG *Img1;
	unsigned int i,j; 
    NumPixNegros = 0;

	Img1=gcGetImgBmp("C01.bmp"); //Lee una imagen en escala de grises

	height 	= Img1->alto;
	width 	= Img1->ancho;

	//Obtiene sus dimensiones
	// printf("\n El alto es: %i",height);
	// printf("\n El ancho es: %i",width);
	// printf("\n El tama~no es: %i",Img1->size);
	printf("\n\n");
	
	//Recorrer los datos de la imagen por completo
	for(i=0; i<Img1->alto; i++) 
	for(j=0; j<Img1->ancho; j++)
		if(Img1->imx[i*Img1->ancho+j] == 0)
			NumPixNegros++;
	printf("Num de pixeles en negro: %i \n", NumPixNegros);

	//* Reserva la memoria para las coordenadas en x & y
	c_x = (int*)malloc(NumPixNegros*sizeof(int));
	if(c_x==NULL){
		printf("Error al reservar memoria");
		exit(0);
	}

	c_y = (int*)malloc(NumPixNegros*sizeof(int));
	if(c_y==NULL){
		printf("Error al reservar memoria");
		exit(0);
	}

	int k=0;
	for(i=0; i<Img1->alto; i++) 
	for(j=0; j<Img1->ancho; j++){
		if(Img1->imx[i*Img1->ancho+j]==0){
			c_x[k]=j;
			c_y[k]=i;
			k++;
		}
	}
	gcFreeImg(Img1);	//* Libera la memoria de la imagen

	int bitsN = contarBits(NumPixNegros); //* Esto es para determinar el numero de bits que se utilizaran 
										//* en los parametros en binario
	printf("El numero %d tiene %d bits.\n", NumPixNegros, bitsN);
}

gcIMG* gcGetImgBmp(char *ruta)
{ gcIMG *img;
  FILE *file;
  int  i,j,a,ar;

// Abrir Archivo de entrada
  if  ( (file = fopen(ruta,"rb"))==NULL )
      { printf(" Error al Abrir Archivo \n");
	exit(1);
	}
// Asignar memoria para la estructura gcIMG
  if  ( (img = (gcIMG *) calloc(1,sizeof(gcIMG)) ) == NULL)
      { printf("Error al reservar memoria para gcIMG \n");
        exit (1);
        }

  fread(img->id,2,1,file);      // Lee 2 bytes del identificador
  fseek(file,10,SEEK_SET);      // Se posiciona en Data offset
  fread(&img->offset,2,1,file); // Lee el offset de la Imagen
  fseek(file,18,SEEK_SET);      // Se posiciona en Width
  fread(&img->ancho,2,1,file);  // Lee el ancho de la Imagen
  fseek(file,22,SEEK_SET);      // Se posiciona en Height
  fread(&img->alto,2,1,file);   // Lee el alto de la Imagen
  fseek(file,28,SEEK_SET);      // Se posiciona en Bits p/pixel
  fread(&img->bpp,1,1,file);    // Lee los Bpp
  fseek(file,34,SEEK_SET);      // Se posiciona en Size
  fread(&img->size,4,1,file);   // Lee el tamaño de la Imagen */

// Comprobar archivo valido
  if  ( (img->id[0]!='B')||(img->id[1]!='M') )
      { printf("Archivo de Formato No Valido \n");
        exit (1);
        }

// Asignar memoria para el encabezado
  if ( (img->head = (unsigned char *) malloc(img->offset)) == NULL )
     { printf("Error al reservar memoria para el encabezado \n");
       exit (1);
       }

// Asignar memoria para la imagen real
  if ( (img->imx =(float *)calloc(img->ancho*img->alto,sizeof(float))) == NULL )
     { printf("Error al reservar memoria para la imagen \n");
       exit (1);
       }

// Lectura del encabezado
  rewind(file);
  fread(img->head,1078,1,file);

// Lectura de la imagen
  a=img->ancho;
  ar=img->size/img->alto;               //calcula el ancho real
  fseek(file,img->offset,SEEK_SET);     // Se posiciona al inicio de la imagen
  for (i=0; i<img->alto; i++)
    { for(j=0; j<img->ancho; j++)
      img->imx[i*a+j]=(float)fgetc(file);
      if(ar!=a) for(j=0;j<ar-a;j++) fgetc(file);  // Si el ancho es mayor
      }                                           // brinca esos datos
  fclose(file);
  img->size=img->ancho*img->alto;       //Asigna el Tamaño Real de la Imagen
  return img;
}

void gcPutImgBmp(char *ruta, gcIMG *img)
{ FILE *file;
  int aux,zero=0,i,j,offset,Newancho;

// Crear un Archivo nuevo
  if ((file = fopen(ruta,"w+b")) == NULL)
        { printf("\nError abriendo el archivo \n");
          exit(1);
          }
//Checar si el ancho es multiplo de 4
  offset=img->ancho%4;
  if (offset) Newancho=img->ancho+(4-offset); //Si no hacerlo multiplo
     else     Newancho=img->ancho;           // Si, si mantenerlo

// Checar el encabezado
  if (img->head) { img->size=(Newancho*img->alto); //Modificar el bitmap size
                   fwrite(img->head,1078,1,file);
                   }
// Generar encabezado:
     else {
            fputc('B',file); fputc('M',file);   // Escribe BMP Identificador
	    aux = Newancho * img->alto + 1078;
	    fwrite(&aux,4,1,file);              // Escribe File Size
	    fwrite(&zero,4,1,file);             // Escribe Word Reserved
            aux=1078;
	    fwrite(&aux,4,1,file);              // Escribe Data Offset
	    // Image Header
	    aux=40;
	    fwrite(&aux,4,1,file);              // Escribe Header Size
	    aux=img->ancho;
	    fwrite(&aux,4,1,file);              // Escribe Width
	    aux=img->alto;
	    fwrite(&aux,4,1,file);              // Escribe Height
	    aux=1;
	    fwrite(&aux,2,1,file);              // Escribe Planes
	    aux=8;
	    fwrite(&aux,2,1,file);              // Escribe Bits p/pixel
	    aux=0;
	    fwrite(&aux,4,1,file);              // Escribe Compression
            aux=(Newancho*img->alto);
	    fwrite(&aux,4,1,file);              // Escribe Bitmap Size
	    aux=0;
	    fwrite(&aux,4,1,file);              // Escribe HResolution
            fwrite(&aux,4,1,file);              // Escribe VResolution
            aux=256;
            fwrite(&aux,4,1,file);              // Escirbe Colors used
            aux=0;
            fwrite(&aux,4,1,file);              // Escirbe Important Colors

// Escritura de la paleta
   	    for (aux=0; aux<256; aux++)
                { for (i=0; i<3; i++) fwrite(&aux,1,1,file);
		  fwrite(&zero,1,1,file);
		  }
	  }
// Escritura del mapa de bits
  aux=img->ancho;
  for(i=0;i<img->alto;i++)
      for(j=0;j<Newancho;j++)
        { if(j>aux-1) fputc(0,file);
          else fputc((unsigned char)img->imx[i*aux+j],file);
         }
  fclose(file);
}

gcIMG *gcNewImg(int ancho,int alto){
  gcIMG *img;
  int i;

  if (( img = (gcIMG *) calloc(1,sizeof(gcIMG)) ) == NULL)
        { printf("Error al reservar memoria para gcIMG\n");
          exit (1);
          }
  img->ancho = ancho;
  img->alto  = alto;
  img->size  = ancho*alto;
  if (( img->imx = (float *) calloc(img->size,sizeof(float)) ) == NULL)
        { printf("Error al reservar memoria para la Imagen \n");
          exit (1);
          }
  img->head = NULL;
  return img;
}

 void gcFreeImg (gcIMG *img)
 { free(img->head);
   free(img->imx);
   free(img);
  }

// void GA::centro_radio(){
//     int A[2], B[2], C[2];

//     for(unsigned int k = 0; k < PobSize; k++){ //* Para toda la poblacion 
//         A[0] = Coorx[Pob[k].Vent[0]];
//         A[1] = Coory[Pob[k].Vent[0]];

//         B[0] = Coorx[Pob[k].Vent[1]];
//         B[1] = Coory[Pob[k].Vent[1]];

//         C[0] = Coorx[Pob[k].Vent[2]];
//         C[1] = Coory[Pob[k].Vent[2]]; 

//         //* Evaluamos que exista el circulo
//         if((A[0] == B[0] && B[0] == C[0]) || (A[1] == B[1] && B[1] == C[1])){
//             Pob[k].cx = 0;
//             Pob[k].cy = 0;
//             Pob[k].r  = 0;
//             return; 
//         }

//         //* Transformacion T
//         double den = ((B[0] - A[0])*(C[1] - A[1])) - ((C[0]-A[0])*(B[1]-A[1])); 
//         den *= 4;
//         // printf("omg: %f", den);

//         double detX =  2*B[0]* pow(C[0], 2) + 2*pow(B[0], 2)*C[1] + 2*pow(B[1], 2)*C[1] - 
//                     2*pow(B[1], 2)*A[1] - 2*pow(A[0], 2)*C[1] - 2*pow(A[1], 2)*C[1]  -
//                     2*pow(C[0], 2)*B[1] + 2*pow(C[0], 2)*A[1] - 2*pow(C[1], 2)*B[1] + 
//                     2*pow(C[1], 2)*A[1] + 2*pow(A[0], 2)*B[1] + 2*pow(A[1], 2)*B[1];
//         detX /= den;
//         // printf("\ndetX: %i", detX);

//         double detY =  2*B[0]*pow(C[0], 2) + 2*B[0]*pow(C[1], 2) - 2*B[0]*pow(A[0], 2) - 
//                     2*B[0]*pow(A[1], 2) - 2*A[0]*pow(C[0], 2) -  2*A[0]*pow(C[1], 2) - 
//                     2*C[0]*pow(B[0], 2) - 2*C[0]*pow(B[1], 2) + 2*C[0]*pow(A[0], 2) + 
//                     2*C[0]*pow(A[1], 2) + 2*A[0]*pow(B[0], 2) + 2*A[0]*pow(B[1], 2);
//         detY /= den;
//         // printf("\ndetY: %i\n", detY);

//         int radio2 = pow((A[0] - detX), 2) + pow((A[1] - detY), 2);
//         int radio = sqrt(radio2);

//         Pob[k].cx = detX;
//         Pob[k].cy = detY;
//         Pob[k].r  = radio;

//         // printf("\ncx: %f", Pob[k].cx);
//         // printf("\ncy: %f", Pob[k].cy);
//         // printf("\nr: %f", Pob[k].r );

//         char a;
//         // std::cin>>a;
//     }
// }


void GA::calcular_centro_radio(int k) {
    double A, B, C, D;
    int Xi, Yi, Xj, Yj, Xk, Yk;

    Xi = Coorx[Pob[k].Vent[0]];
    Yi = Coory[Pob[k].Vent[0]];
    Xj = Coorx[Pob[k].Vent[1]];
    Yj = Coory[Pob[k].Vent[1]];
    Xk = Coorx[Pob[k].Vent[2]];
    Yk = Coory[Pob[k].Vent[2]];

    if ((Xi == Xj && Xi == Xk) || (Yi == Yj && Yj == Yk)) {
        Pob[k].cx = 0;
        Pob[k].cy = 0;
        Pob[k].r = 0;
    } else {
        calcular_coeficientes(A, B, C, D, Xi, Yi, Xj, Yj, Xk, Yk, Pob[k].cx, Pob[k].cy);
        Pob[k].r = calcular_distancia(Pob[k].cx, Pob[k].cy, Coorx[Pob[k].Vent[0]], Coory[Pob[k].Vent[0]]);
    }
}

void GA::calcular_coeficientes(double &A, double &B, double &C, double &D, int Xi, int Yi, int Xj, int Yj, int Xk, int Yk, double &cx, double &cy) {
    A = pow(Xj, 2) + pow(Yj, 2) - (pow(Xi, 2) + pow(Yi, 2));
    B = 2 * (Yj - Yi);
    C = pow(Xk, 2) + pow(Yk, 2) - (pow(Xi, 2) + pow(Yi, 2));
    D = 2 * (Yk - Yi);

    cx = ((A * D) - (B * C)) / (4 * ((Xj - Xi) * (Yk - Yi) - (Xk - Xi) * (Yj - Yi)));

    A = 2 * (Xj - Xi);
    B = pow(Xj, 2) + pow(Yj, 2) - (pow(Xi, 2) + pow(Yi, 2));
    C = 2 * (Xk - Xi);
    D = pow(Xk, 2) + pow(Yk, 2) - (pow(Xi, 2) + pow(Yi, 2));

    cy = ((A * D) - (B * C)) / (4 * ((Xj - Xi) * (Yk - Yi) - (Xk - Xi) * (Yj - Yi)));
}

double GA::calcular_distancia(double x1, double y1, double x2, double y2) {
    return sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
}

void GA::centro_radio() {
    for (int k = 0; k < PobSize; k++) {
        calcular_centro_radio(k);
    }
}


//? Parametros del algoritmo genetico (Segun el archivo presentado)
const unsigned int  NumeroDeGeneraciones = 3000;
const unsigned int  SizePoblacion        = 70;
const unsigned int  NumeroDeGenes        = 3; //* x, y, r

const double        ProbCruza = 0.55;
const double        ProbMuta  = 0.10;

int main(){
    srand(time(0));
	read();
	unsigned int nb = contarBits(NumPixNegros);

	const unsigned int  NumeroDeBitsxGen[NumeroDeGenes] = {nb, nb, nb};    

    GA ga(NumeroDeGeneraciones, SizePoblacion, NumeroDeGenes, NumeroDeBitsxGen, NumPixNegros, c_x, c_y, height, width);
    printf("\nSe genero el GA...\n");
    
    ga.process(ProbCruza, ProbMuta);

    return 0;
}

