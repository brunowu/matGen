#include "edmg.h"
#include "ebmg.h"

static char help[] = "Dense Matrix Generator by select eigenvalues";

int main(int argc, char ** argv){
	 int             world_size;
	 PetscBool       edmg_flg;

	PetscInitialize(&argc,&argv,(char *)0,help);

	MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
    	PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "Using number of %d precessors for the generation...\n", world_size);
	PetscOptionsGetBool(NULL, PETSC_NULL, "-dense", &edmg_flg, NULL);

	if(edmg_flg){
		EDMG();
	} 
	else {
		EBMG();
	}

	PetscFinalize();

	return 0;
}
