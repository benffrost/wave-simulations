//CODE SAMPLE FOR INCLUSION WITH BEN FROST'S RESUME
__global__ void TimeStep_DtoE_TE_Bloch(			cuDoubleComplex*  Ex, cuDoubleComplex*  Ey,
												const cuDoubleComplex*  Dx, const cuDoubleComplex*  Dy,
												const matrix22* Eps,
												const int2 dims, const cuDoubleComplex blochfactorx, const cuDoubleComplex blochfactory)
/*
anisotropic epsilon materials handling for FDTD on the Yee grid based on the "WCMOD7" algorithm described in:
http://www.sciencedirect.com/science/article/pii/S0021999107002227
Scrapes one dimension of threads over the two-dimensional arrays; expects a four-pixel buffer around the fields.
Takes the D field and applies the inverted epsilon tensor to make the E field.
if needed, we can use this for this for BtoH in TE because the operation works the same even if the letters aren't (but not for non-bloch cases)

note that we don't use vectors to store the fields, as that would be slower :( and is potentially even more confusing to read; in particular it would be inefficient/unruly for non-periodic cases, although this isn't a big deal here.
also note that this algorithm may feature errors if the cuda blocks and simulation boundaries line up exactly, but fixing this would slow the algorithm down more than simply requiring that the boundaries don't line up perfectly.

apologies for the complex math handling's lack of operators, but this is written to use the default cuComplex.h library
*/
{
	const int gtidx 	= blockDim.x * blockIdx.x + threadIdx.x;
	const int workx_dx 	= threadIdx.x + 1;
	const int workx_dy 	= threadIdx.x + 0;
	int index 			= gtidx + 2 * dims.x;

	const bool validw = ((gtidx < dims.x - 2) && (gtidx > 1));	 
	
	__shared__  cuDoubleComplex Dxtile[blockwidth_max + 1];
	__shared__  cuDoubleComplex Euytile[blockwidth_max + 1];

	cuDoubleComplex Eux_Above      = COMPLEXZERO;
	cuDoubleComplex Eux_Here       = COMPLEXZERO;
	
	cuDoubleComplex Dy_Here        = COMPLEXZERO;
	cuDoubleComplex Dy_Below       = COMPLEXZERO;

	cuDoubleComplex Dy_Here_halo   = COMPLEXZERO;
	cuDoubleComplex Dy_Below_halo  = COMPLEXZERO;

///////////////presteps:
//DY PRESTEP
if(validw) //indexing all of this is crazy
	{
		if(threadIdx.x == 0)
			Dxtile[workx_dx - 1]		= Dx[index - 1];	
		else if(threadIdx.x == blockDim.x- 1)	
			Dy_Here_halo 				= cuCmul(Dy[index + 1 + dims.x * (dims.y - 5)],cuCconj(blochfactory));						

		Dy_Here 						= cuCmul(Dy[index + dims.x * (dims.y - 5)],cuCconj(blochfactory));;
		Dxtile[workx_dx]				= Dx[index];
	}
	else if(gtidx == 1)
		Dxtile[workx_dx]				= cuCmul(Dx[index + dims.x - 4],cuCconj(blochfactorx)); 
	else if(gtidx == dims.x - 2)
		Dy_Here 						= cuCmul(cuCmul(Dy[index + dims.x * (dims.y - 6) + 4],cuCconj(blochfactory)),blochfactorx);

	__syncthreads();
	if(validw)
		Eux_Above	= cuCscale(cuCadd(Dxtile[workx_dx - 1],Dxtile[workx_dx]), 0.5 * Eps[index].yx); 	
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////END PRELOADS
#pragma UNROLL 9
	for(int iy = 2; iy < dims.y - 3; iy++)
	{
		Eux_Here 		= Eux_Above;
		Dy_Below       	= Dy_Here;
			
		if(validw)
		{
			if(threadIdx.x == 0)
				Dxtile[workx_dx - 1]	= Dx[index - 1 + dims.x];
			else if(threadIdx.x == blockDim.x- 1)	
			{
				Dy_Below_halo 			= Dy_Here_halo;
				Dy_Here_halo 			= Dy[index + 1];
				Euytile[workx_dy + 1] 	= cuCscale(cuCadd(Dy_Here_halo,Dy_Below_halo), 0.5 * Eps[index + 1].xy);					
			}
			
			Dy_Here 				= Dy[index];
			Euytile[workx_dy] 		= cuCscale(cuCadd(Dy_Here,Dy_Below), 0.5 * Eps[index].xy);	
			Dxtile[workx_dx]		= Dx[index + dims.x];
		}
		else if(gtidx == 1)
			Dxtile[workx_dx]				= cuCmul(Dx[index + 2 * dims.x - 4],cuCconj(blochfactorx)); 
		else if(gtidx == dims.x - 2)
		{
			Dy_Here 					= cuCmul(Dy[index + 4 - dims.x],blochfactorx);
			Euytile[workx_dy] 			= cuCscale(cuCadd(Dy_Here,Dy_Below), 0.5 * Eps[index].xy);	
		}
		///
		__syncthreads();
		if(validw)
		{
			Eux_Above	= cuCscale(cuCadd(Dxtile[workx_dx - 1],Dxtile[workx_dx]), 0.5 * Eps[index + 1 * dims.x].yx); 	
			Ex[index] 	= cuCadd(cuCscale(Dx[index],0.5 * (Eps[index].xx + Eps[index + 1].xx)),cuCscale(cuCadd(Euytile[workx_dy],Euytile[workx_dy + 1]),0.5));
			Ey[index] 	= cuCadd(cuCscale(Dy_Here,0.5 * (Eps[index].yy + Eps[index + dims.x].yy)),cuCscale(cuCadd(Eux_Above,Eux_Here),0.5));	
		}
		index += dims.x;
		__syncthreads();
	}

	/////////////////////////wrapup
	Eux_Here 	= Eux_Above;
	Dy_Below    = Dy_Here;
		
	if(validw)
	{
		if(threadIdx.x == 0)
			Dxtile[workx_dx - 1]	= cuCmul(Dx[index - 1 + dims.x * (5 - dims.y)],blochfactory);	
		else if(threadIdx.x == blockDim.x- 1)
		{
			Dy_Below_halo 			= Dy_Here_halo;
			Dy_Here_halo 			= Dy[index + 1];
			
			Euytile[workx_dy + 1] 	= cuCscale(cuCadd(Dy_Here_halo,Dy_Below_halo), 0.5 * Eps[index + 1].xy);					
		}		
		Dy_Here 					= Dy[index];
		Euytile[workx_dy] 			= cuCscale(cuCadd(Dy_Here,Dy_Below), 0.5 * Eps[index].xy);	

		Dxtile[workx_dx]			= cuCmul(Dx[index + dims.x * (5 - dims.y)],blochfactory);
	}
	else if(gtidx == 1)
		Dxtile[workx_dx]			=  cuCmul(cuCmul(Dx[index - 4 + dims.x * (6 - dims.y)],blochfactory), cuCconj(blochfactorx)); 
	else if(gtidx == dims.x - 2)
	{
		Dy_Here 					=  cuCmul(Dy[index + 4 - dims.x],blochfactorx);
		Euytile[workx_dy] 			= cuCscale(cuCadd(Dy_Here,Dy_Below), 0.5 * Eps[index].xy);	
	}

	///
	__syncthreads();
	if(validw)
	{
		Eux_Above	= cuCscale(cuCadd(Dxtile[workx_dx - 1],Dxtile[workx_dx]), 0.5 * Eps[index].yx); 	
		Ex[index] 	= cuCadd(cuCscale(Dx[index],Eps[index].xx),cuCscale(cuCadd(Euytile[workx_dy],Euytile[workx_dy + 1]),0.5));
		Ey[index] 	= cuCadd(cuCscale(Dy_Here,Eps[index].yy),cuCscale(cuCadd(Eux_Above,Eux_Here),0.5));	
	}	
}