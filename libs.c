#include "libs.h"

PetscErrorCode getFileSize(const char * name, PetscInt * size){
  FILE * fptr;
  *size = 0L;

#ifdef LINUX
  struct stat fs;

  if(stat(name,&fs)!=0){
    perror("Cannot state file\n");
  }
  *size=fs.st_size;

#else
fptr=fopen(name,"rb");
  if(fptr!=NULL){
    fseek(fptr,0L,SEEK_END);
    *size = ftell(fptr);
    fclose(fptr);
  }
#endif

  return 0;
}

PetscErrorCode readBinaryScalarArray(const char * name, PetscInt * nb, PetscScalar * array){
  int file_descriptor;
  PetscErrorCode ierr;
  PetscInt size;

  getFileSize(name,&size);

  if(*nb<=0) *nb=(PetscInt)size/((PetscInt)sizeof(PetscScalar));
  if(size/sizeof(PetscScalar)!=*nb) {
    return 1;
  }


  ierr=PetscBinaryOpen(name,FILE_MODE_READ,&file_descriptor);CHKERRQ(ierr);
  ierr=PetscBinarySynchronizedRead(PETSC_COMM_WORLD,file_descriptor,array,*nb,PETSC_SCALAR);CHKERRQ(ierr);
  ierr=PetscBinaryClose(file_descriptor);CHKERRQ(ierr);

  return ierr;
}


void random_selection(PetscScalar *ret, PetscInt nombre)
{

  PetscInt	i;

  int 			my_seed,my_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&my_rank);
  my_seed=time(NULL)+my_rank;

  srand(my_seed); 

  for(i = 0; i < nombre; i++){

/*
    if(i<1*nombre/4) ret = 1*(cos(2*PI*i/nombre)+2) + PETSC_i*1*(sin(2*PI*i/nombre)+2);
    else if(i<2*nombre/4) ret =1* (cos(2*PI*i/nombre)+2) + PETSC_i*1*(sin(2*PI*i/nombre)-2);
    else if(i<3*nombre/4) ret = 1*(cos(2*PI*i/nombre)+4) + PETSC_i*1*(sin(2*PI*i/nombre)+2);
    else ret = 1*(cos(2*PI*i/nombre)+4) + PETSC_i*1*(sin(2*PI*i/nombre)-2);
*/

/*
    if(i < 9*nombre / 10)
    	ret[i] = -1*(cos(2*PI*i/nombre)+2+5*rand()/RAND_MAX) + PETSC_i*1*(sin(2*PI*i/nombre));
    else ret[i] =-1*(cos(2*PI*i/nombre)+2) + PETSC_i*1*(sin(2*PI*i/nombre) );
*/

/*

    if(i < 1*nombre / 10)
        ret[i] = 1*(cos(2*PI*i/nombre)+1.2500) + PETSC_i*2*(sin(2*PI*i/nombre));
    else ret[i] =1*(cos(2*PI*i/nombre)+1.250) + PETSC_i*2*(sin(2*PI*i/nombre));

*/

    if(i < 5*nombre / 10)
        ret[i] = 5*(cos(PI*i/nombre)) - 4.2  + PETSC_i*3*(sin(PI*i/nombre));
    else ret[i] =5*(cos(PI*i/nombre)) + 4.2  + PETSC_i*3*(sin(PI*i/nombre));


/*
    if(i < 3*nombre / 10)
        ret[i] = 1*(cos(2*PI*i/nombre)+2) + PETSC_i*1*sin(2*PI*i/nombre);
    else if(i < 6 * nombre / 10) ret[i] = -10*(cos(2*PI*i/nombre)+200) + PETSC_i*1*(sin(2*PI*i/nombre) + 277);
    else ret[i] = 1*(cos(2*PI*i/nombre)+3000) + PETSC_i*1*(sin(2*PI*i/nombre) - 2800); 
*/
    }

}

void selection(PetscScalar *ret, PetscInt nombre, PetscInt min, PetscInt max)
{

  PetscScalar   *tab;
  PetscInt    i, indice, maxi = max - min;
 
  if(min >= max || nombre > maxi + 1 || nombre < 1)
    PetscPrintf(PETSC_COMM_WORLD,"Input values of tirage() are wrong.\n");
 
  PetscMalloc((maxi + 1)*sizeof(tab[0]),&tab);

  for(i = 0; i < maxi + 1; i++)
    tab[i] = i + min;
 
  for(i = 0; i < nombre; i++){
    indice = rand() % (maxi + 1);
    ret[i] = tab[indice];
    tab[indice] = tab[maxi];
    maxi--;
  }
  PetscFree(tab);
}


void change(PetscScalar *array, PetscInt n, PetscReal ratio)
{
        PetscInt i;

        for(i = n-1; i>0;i--){
                if(i < ratio * n){
			array[i] = -array[i];
		}
                else{
			array[i] = array[i];
		}
        }
}

void shuffer(PetscScalar *array, PetscInt n)
{
	PetscInt index, i;
	PetscScalar tmp;
	srand(time(NULL));

	for(i = n-1; i>0;i--){
		index = rand() % i;
		tmp = array[i];
		array[i]=array[index];
		array[index]=tmp;
	}
}

PetscInt *indexShuffer(PetscInt n)
{
	PetscInt index, i;
	PetscInt *a;
	a = (PetscInt *) malloc(n*sizeof(PetscInt));
	PetscInt tmp;
	srand(time(NULL));
	for (i = 0; i < n; i++)
		a[i] = i;

	for(i = n-1; i>0;i--){
		index = rand() % i;
		tmp = a[i];
		a[i]=a[index];
		a[index]=tmp;
	}
	return a;
}

PetscInt factorial(PetscInt start, PetscInt end) { 
  PetscInt i;
  PetscInt valeur;
  if(start>end){
    valeur = 1;
  }else{
    valeur = start;
    for(i= start+1; i <= end; i++){
        valeur *= i;
    }
  }
  return valeur;
}

PetscScalar Random (PetscInt _iMin, PetscInt _iMax) 
{ 
  return (_iMin + (rand () % (_iMax-_iMin+1))); 
} 

PetscInt IRandom (PetscInt _iMin, PetscInt _iMax) 
{ 
  return (_iMin + (rand () % (_iMax-_iMin+1))); 
} 

void printarray(PetscInt n, PetscScalar *a) {

	PetscScalarView(n, a, PETSC_VIEWER_STDOUT_WORLD);
}
