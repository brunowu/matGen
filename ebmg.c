#include "ebmg.h"

PetscErrorCode EBMG(){

	char fileb[PETSC_MAX_PATH_LEN];

	Vec            eigenvalues;
  	Mat            Mt, AM, MA, matAop;
  	Mat            A;
  	PetscErrorCode ierr;
  	PetscInt       n,i,j,k,degree;
  	PetscScalar    rRandom1, rRandom2;
  	PetscInt       iRandom, d1, d2;
  	PetscInt    	size;
  	PetscBool      flagb,n_flg,degree_flg,d1_flg,d2_flg;
	MatInfo     	Ainfo;
	double        	gnnz;
	char           	matrixOutputFile[PETSC_MAX_PATH_LEN];
	PetscViewer    	output_viewer;
	clock_t start, finish;
	double  duration;

	PetscInt Istart, Iend;

    ierr=PetscOptionsGetInt(NULL, PETSC_NULL,"-n",&n, &n_flg);CHKERRQ(ierr);
	ierr=PetscOptionsGetInt(NULL, PETSC_NULL,"-degree",&degree, &degree_flg);CHKERRQ(ierr);
	ierr=PetscOptionsGetInt(NULL, PETSC_NULL,"-d1",&d1, &d1_flg);CHKERRQ(ierr);
	ierr=PetscOptionsGetInt(NULL, PETSC_NULL,"-d2",&d2, &d2_flg);CHKERRQ(ierr);

    if (!n_flg){
		ierr = PetscPrintf(PETSC_COMM_WORLD,"!!!Please set the dimension of matrix to be generated\n");CHKERRQ(ierr);
		return 0;
	}

	if (!degree_flg){
		degree = 5;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Unset degree, using default degree = %d\n", degree);CHKERRQ(ierr);
	}
	
	if (!d1_flg){
		d1 = 0;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Unset d1, using default degree = %d\n", d1);CHKERRQ(ierr);
	}

	if (!d2_flg){
		d2 = 50;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Unset d2, using default degree = %d\n", d2);CHKERRQ(ierr);
	}

    PetscPrintf(PETSC_COMM_WORLD, "To generate matrix with dim = %d, degree = %d, d1 = %d, d2 = %d ... \n", n, degree,d1,d2);

	ierr = VecCreate(PETSC_COMM_WORLD,&eigenvalues);CHKERRQ(ierr); 
	ierr = VecSetSizes(eigenvalues,PETSC_DECIDE,n);CHKERRQ(ierr); 
	ierr = VecSetFromOptions(eigenvalues);CHKERRQ(ierr);

	PetscInt istart, iend;

	ierr = VecGetOwnershipRange(eigenvalues, &istart, &iend);

	PetscScalar *Deigenvalues;
  	PetscMalloc1(n,&Deigenvalues);

 	PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");
    ierr=PetscOptionsGetString(NULL,PETSC_NULL,"-vfile",fileb,PETSC_MAX_PATH_LEN-1,&flagb);CHKERRQ(ierr);

	if (!flagb){
		PetscPrintf(PETSC_COMM_WORLD, "Not providing the outside eigenvalues files, using the internal functions to generate them...\n");
		random_selection(Deigenvalues,n);
		shuffer(Deigenvalues,n);
	}
	else{
        PetscPrintf(PETSC_COMM_WORLD, "Using the eigenvalues provides by outside files: %s ...\n", fileb);
		readBinaryScalarArray(fileb, &size, Deigenvalues);
		change(Deigenvalues, size, 0.5);
		shuffer(Deigenvalues,size);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"read file size = %ld\n", size);CHKERRQ(ierr);
		if(size != n){
			ierr = PetscPrintf(PETSC_COMM_WORLD,"!!!read file size and vec dimemson do not match and mat dim set to be equal to vec dim\n");CHKERRQ(ierr);
			return 0;
		}
		n = size;
	}
	PetscScalar tmp;

	for(i=istart;i<iend;i++){
		tmp = Deigenvalues[i];
		VecSetValue(eigenvalues, i, tmp, INSERT_VALUES);

	}

	start = clock();
	ierr = MatCreate(PETSC_COMM_WORLD,&Mt);CHKERRQ(ierr);
	ierr = MatSetSizes(Mt,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
	ierr = MatSetType(Mt,MATMPIAIJ);CHKERRQ(ierr);
	ierr = MatSetFromOptions(Mt);CHKERRQ(ierr); 
	ierr = MatSetUp(Mt);CHKERRQ(ierr);
	ierr = MatGetOwnershipRange(Mt,&Istart,&Iend); CHKERRQ(ierr);
	ierr = MatDiagonalSet(Mt,eigenvalues,INSERT_VALUES);CHKERRQ(ierr);

	for (i=Istart; i<Iend; i++){
	  for (j=i-d1; j<i; j++){
	  	if(j >= 0)
	  	{
	    	rRandom1 = Random (0,10);
	    	ierr = MatSetValues(Mt,1,&i,1,&j,&rRandom1,INSERT_VALUES);CHKERRQ(ierr);
	  	}
	   }
	}

	for (i=Istart; i<Iend; i++){
	  for (j=i; j<i+d2; j++){
	  	if(j < n){
	    	rRandom2 = Random (0,10);
	    	ierr = MatSetValues(Mt,1,&i,1,&j,&rRandom2,INSERT_VALUES);CHKERRQ(ierr);
		}
	   }
	}

	ierr = MatAssemblyBegin(Mt,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Mt,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatDuplicate(Mt,MAT_DO_NOT_COPY_VALUES,&A);CHKERRQ(ierr);

	PetscInt iI, jJ;
	PetscScalar vV=1;
	i=0;

	for (i=Istart; i<Iend; i++){
		if (i<degree) {
	  		iI=i;
	  		jJ=i+1;
	  		ierr = MatSetValues(A,1,&iI,1,&jJ,&vV,INSERT_VALUES);CHKERRQ(ierr);
	  		i++;
		}
	}	

	i=degree+1;

	for (i=Istart; i<Iend; i++){
		if (i<n-1) {
	  		srand (time (NULL)); 
	  		iRandom = IRandom (1,degree);
	  		for(j=0; j<min(iRandom,n-1-i); j++)
	  		{
	    		iI=i+j;
	    		jJ=i+j+1;
	    		ierr = MatSetValues(A,1,&iI,1,&jJ,&vV,INSERT_VALUES);CHKERRQ(ierr);
	  		}

	  		i=i+min(iRandom,n-1-i);
	  		srand (time (NULL)); 
	  		iRandom = IRandom (1,max(1,n-1-i));
	  		i=i+iRandom;
		}
	}

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
	PetscInt my_factorielle_bornes=1;
	ierr = MatDuplicate(Mt,MAT_COPY_VALUES,&matAop);CHKERRQ(ierr);

	my_factorielle_bornes =  factorial(1,(2*degree-2));
	ierr = MatScale(Mt,my_factorielle_bornes);CHKERRQ(ierr);

	for (k=1; k<=(2*degree-2); k++) {
	  ierr = MatMatMultSymbolic(matAop,A,PETSC_DEFAULT,&MA);CHKERRQ(ierr);
	  ierr = MatMatMultNumeric(matAop,A,MA);CHKERRQ(ierr);

	  ierr = MatMatMultSymbolic(A,matAop,PETSC_DEFAULT,&AM);CHKERRQ(ierr);
	  ierr = MatMatMultNumeric(A,matAop,AM);CHKERRQ(ierr);

	  ierr = MatAYPX(matAop,0,AM,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	  ierr = MatAXPY(matAop,-1,MA,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	  my_factorielle_bornes =  factorial(k+1,2*degree-2);
	  ierr = MatAXPY(Mt,my_factorielle_bornes,matAop,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	  ierr = MatZeroEntries(AM);CHKERRQ(ierr);
	  ierr = MatZeroEntries(MA);CHKERRQ(ierr);
	}
	
	MatGetInfo(Mt,MAT_GLOBAL_SUM,&Ainfo);
	gnnz = Ainfo.nz_used;

	sprintf(matrixOutputFile,"EBMG_matrix_nb_%d_%dx%d_%g_nnz.gz",n, n,n,gnnz);
		
	PetscPrintf(PETSC_COMM_WORLD,"\n@>Dumping matrix to PETSc binary %s\n",matrixOutputFile);
			
	PetscViewerBinaryOpen(PETSC_COMM_WORLD,matrixOutputFile,FILE_MODE_WRITE,&output_viewer);
	PetscViewerPushFormat(output_viewer,PETSC_VIEWER_ASCII_INFO_DETAIL);
	//MatView(Mt,PETSC_VIEWER_STDOUT_WORLD);
	MatView(Mt,output_viewer);
	PetscViewerDestroy(&output_viewer);
		
	PetscPrintf(PETSC_COMM_WORLD,"\n@>Matrix %s Dumped\n\n",matrixOutputFile);
	PetscPrintf(PETSC_COMM_WORLD,"\n>>>>>>Please use the command 'gzip -d **' to unzip the file to binary file\n\n");

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
    PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD,"\nElapsed time is %f seconds\n\n", duration);
    PetscPrintf(PETSC_COMM_WORLD, "--------------------------\n");

	return ierr;

}
