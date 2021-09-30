#include "DLLFunctions_new.cuh"

int main(int argc, char *argv[])
{
	if (argc > 1)
	{
		FILE* ConfigFile = fopen(argv[1], "rt");

		long numX, numY, numx, numy, numF;
		
		long Cent1X, Cent1Y, Cent2X, Cent2Y;
		float Rad1, Rad2;
		uint8_t DFEh = 0;

		char BFileName[512], DataFileName[512], OutputFileName_P1[512], OutputFileName_P2[512], OutputFileName_Raw[512];

		fscanf(ConfigFile, "%s", &DataFileName[0]);
		fscanf(ConfigFile, "%s", &BFileName[0]);
		fscanf(ConfigFile, "%s", &OutputFileName_P1[0]);
		fscanf(ConfigFile, "%s", &OutputFileName_P2[0]);
		fscanf(ConfigFile, "%s", &OutputFileName_Raw[0]);
		fscanf(ConfigFile, "%d", &numX);
		fscanf(ConfigFile, "%d", &numY);
		fscanf(ConfigFile, "%d", &numx);
		fscanf(ConfigFile, "%d", &numy);
		fscanf(ConfigFile, "%d", &numF);
		fscanf(ConfigFile, "%d", &Cent1X);
		fscanf(ConfigFile, "%d", &Cent1Y);
		fscanf(ConfigFile, "%f", &Rad1);
		fscanf(ConfigFile, "%d", &Cent2X);
		fscanf(ConfigFile, "%d", &Cent2Y);
		fscanf(ConfigFile, "%f", &Rad2);

		fclose(ConfigFile);

		uint8_t* Background = (uint8_t*)malloc(sizeof(uint8_t) * numX * numY * 3 / 2);

		if (strcmp(BFileName, "INVALID") != 0)
		{
			FILE* BFID = fopen(BFileName, "rb");
			fread(Background, sizeof(uint8_t), numX * numY * 3 / 2, BFID);
			fclose(BFID);
		}
		else
		{
			long idx = 0;
			for (idx = 0; idx < (numX * numY * 3 / 2); idx++)
			{
				Background[idx] = 0u;
			}
		}
		printf("----------------------------------------------------\n  SPOOF OCM POST PROCESSING LOG FILE  \n----------------------------------------------------\n");

		printf("File name:        %s\n", DataFileName);
		printf("Background name:  %s\n", BFileName);
		printf("Output P1 name:   %s\n", OutputFileName_P1);
		printf("Output P2 name:   %s\n", OutputFileName_P2);
		printf("%d x %d --> %d x %d (%d frames)\n", numX, numY, numx, numy, numF);
		printf("(%d,%d, r = %1.1f) and (%d,%d, r = %1.1f)\n", Cent1X, Cent1Y, Rad1, Cent2X, Cent2Y, Rad2);

		cudaError_t Ret1 = FFPSOCM_InitializeEverything(numX, numY, numx, numy, Background, &Rad1, &Rad2);
		printf("Initialized (%s)\n", cudaGetErrorString(Ret1));

		uint8_t* Input = (uint8_t*)malloc(sizeof(uint8_t) * numX * numY * 3 / 2);
		cufftComplex* RawIm = (cufftComplex*)malloc(sizeof(cufftComplex) * numX * numY);
		cufftComplex* Im1 = (cufftComplex*)malloc(sizeof(cufftComplex) * numx * numy);
		cufftComplex* Im2 = (cufftComplex*)malloc(sizeof(cufftComplex) * numx * numy);

		long fidx = 0;
		FILE* DFID = fopen(DataFileName, "rb");
		FILE* P1FID = fopen(OutputFileName_P1, "wb");
		FILE* P2FID = fopen(OutputFileName_P2, "wb");
		FILE* RawFID = NULL;
		if (WRITERAWOUT) RawFID = fopen(OutputFileName_Raw, "wb");
		clock_t t = clock(), t2;
		for (fidx = 0; fidx < numF; fidx++)
		{
			fread(Input, sizeof(uint8_t), (numX * numY * 3) / 2, DFID);

			cudaError_t Ret2 = FFPSOCM_ProcessFrame(Input, RawIm, Im1, Im2, numX, numY, numx, numy, &Cent1X, &Cent1Y, &Cent2X, &Cent2Y,
				&Rad1, &Rad2, &DFEh);

			fwrite(Im1, sizeof(cufftComplex), numx * numy, P1FID);
			if (WRITERAWOUT) fwrite(RawIm, sizeof(cufftComplex), numX * numY, RawFID);
			fwrite(Im2, sizeof(cufftComplex), numx * numy, P2FID);
			t2 = clock() - t;
			printf("\rFinished %d/%d frames in %1.1f seconds...                        ", fidx, numF, ((float)t2) / CLOCKS_PER_SEC);
		}

		cudaError_t Ret3 = FFPSOCM_DestroyEverything();
		printf("\nDestroyed everything (%s)\n", cudaGetErrorString(Ret3));
		fclose(DFID);
		fclose(P1FID);
		fclose(P2FID);
		if (WRITERAWOUT) fclose(RawFID);

		free(Input);
		free(RawIm);
		free(Background);
		free(Im1);
		free(Im2);
	}
	return 0;
}