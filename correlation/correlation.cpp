/*********************************************************\
 * UNREDD+UNREDD+UNR     NREDD+U            +UNREDD      *
 * DD+UNREDD+UNREDD+     D+UNRED         UNREDD+UNREDD   *
 * NREDD+UNREDD+UNRE    NREDD+UNR       EDD+UNREDD+UNRE  *
 * D+UNR                D+UNREDD+      +UNRED     EDD+UN *
 * REDD+               NREDD UNRED     EDD+U       NREDD *
 * +UNRE               D+UNR DD+UN     UNRED       D+UNR *
 * EDD+UNREDD+UNRE    NREDD+UNREDD+    DD+UN       REDD+ *
 * UNREDD+UNREDD+U    D+UNREDD+UNRE    NREDD       +UNRE *
 * DD+UNREDD+UNRED   NREDD+UNREDD+UN   D+UNR       EDD+U *
 * NREDD             D+UNR     NREDD   REDD+       UNRED *
 * D+UNR            NREDD       +UNRE  +UNRED     EDD+UN *
 * REDD+            D+UNR       EDD+U   DD+UNREDD+UNRED  *
 * +UNRE           NREDD         NREDD   REDD+UNREDD+U   *
 * EDD+U           D+UNR         D+UNR      REDD+UN      *
\********** Implemented by Dr. Matthias Schramm **********/
//ver 20120522


/*********************************************************\
 * PROJECT: CORRELATION                                  *
 * Goal:                                                 *
 * - Calculate the cross-correlation of two images       *
 * - Get information about relative georectification     *
 * Input:                                                *
 * - Path and name of metadata file                      *
 * - all needed information in metadata file             *
 * - path and names of input images                      *
 * - needed output format                                *
 * - calculated pixels of input images                   *
 * Output:                                               *
 * - depending on output format                          *
 * - mean x- and y-gap for whole images                  *
 * - x- and y-gap for each pixels of images overlap      *
 * - x- and y-gaps for defined pixels                    *
\*********************************************************/

#include "correlation.h"

class Parallelizing
/*********************************************************\
 * CLASS: PARALLELIZING                                  *
 * Goal:                                                 *
 * - Parallelize process of cross-correlation            *
 * - Calculate several pixels at the same time           *
 * - 100% CPU                                            *
 * Input:                                                *
 * - Arrays of the two images (only used overlap pixels) *
 * - Array of output image                               *
 * - size of overlap images (number of pixels)           *
 * - size of moving window and maximal gap (in pixels)   *
 * - counters for calculated pixels (for statistics)     *
 * Output:                                               *
 * - filled output array                                 *
\*********************************************************/
{
	float ***image1, ***image2, ***image_invar, ***output;
	int bandsDataset1, bandsDataset2;
	int no_pix_x,no_pix_y,window_size,max_gap,value_invar;
	long *counter_abs, *counter_used;
public:
/*********************************************************\
 * FUNCTION: OPERATOR                                    *
 * Goal:                                                 *
 * - Parallelizing                                       *
 * Input:                                                *
 * - range of parallelizing                              *
 * Output: void                                          *
\*********************************************************/
	void operator()(const blocked_range<size_t>& r) const
	{
		float counter=0.0;
		for (size_t runtime=r.begin(); runtime!=r.end(); ++runtime)
		{
			long var1=*counter_abs;
			*counter_abs=var1+1;
			bool break_=0;
			int centre_pix_y= (int)(runtime/no_pix_x)+(int)((window_size+1)/2)+max_gap-1;
			int centre_pix_x= (int)(runtime%no_pix_x)+(int)((window_size+1)/2)+max_gap-1;	

/*********************************************************\
 * Testing, if actual pixel is invariant (else not used) *
\*********************************************************/
			if(image_invar[0][centre_pix_x][centre_pix_y] != value_invar)
			{
				output[0][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
				output[1][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
				output[2][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
				break_=1;
			}
/* if calculation aborted: stop actual window calc.      */
			if (break_==1)
			{
				continue;
			}

/*********************************************************\
 * Testing, if actual pixel window can be used           *
 * if background pixel (0) in window or gap, not usable  *
\*********************************************************/
//--------------------------------------------------------
// Image 1
//--------------------------------------------------------
/* define and declare new moving window and gap Matrix   */
			MatrixXd* actWindowAndGap1 = new MatrixXd[bandsDataset1];
			for (int k=0; k<bandsDataset1; k++)
				actWindowAndGap1[k].resize(window_size+(2*max_gap),window_size+(2*max_gap));
/* fill new matrix with input data                       */
			for (int k=0; k<bandsDataset1; k++)
				for (int i=centre_pix_x-(int)(window_size/2)-max_gap; i<=centre_pix_x+(int)(window_size/2)+max_gap; i++)
					for (int j=centre_pix_y-(int)(window_size/2)-max_gap; j<=centre_pix_y+(int)(window_size/2)+max_gap; j++)
					{
						actWindowAndGap1[k](i-(centre_pix_x-(int)(window_size/2)-max_gap),j-(centre_pix_y-(int)(window_size/2)-max_gap)) = image1[k][i][j];
						if(actWindowAndGap1[k](i-(centre_pix_x-(int)(window_size/2)-max_gap),j-(centre_pix_y-(int)(window_size/2)-max_gap)) == 0)
						{
							output[0][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
							output[1][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
							output[2][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
							break_=1;
							break;
						}
					}
/* if calculation aborted: stop actual window calc.      */
			if (break_==1)
			{
				delete[] actWindowAndGap1;
				continue;
			}

//--------------------------------------------------------
// Image 2
//--------------------------------------------------------
/* define and declare new moving window and gap Matrix   */
			MatrixXd* actWindowAndGap2 = new MatrixXd[bandsDataset2];
			for (int k=0; k<bandsDataset2; k++)
				actWindowAndGap2[k].resize(window_size+(2*max_gap),window_size+(2*max_gap));
/* fill new matrix with input data                       */
			for (int k=0; k<bandsDataset2; k++)
				for (int i=centre_pix_x-(int)(window_size/2)-max_gap; i<=centre_pix_x+(int)(window_size/2)+max_gap; i++)
					for (int j=centre_pix_y-(int)(window_size/2)-max_gap; j<=centre_pix_y+(int)(window_size/2)+max_gap; j++)
					{
						actWindowAndGap2[k](i-(centre_pix_x-(int)(window_size/2)-max_gap),j-(centre_pix_y-(int)(window_size/2)-max_gap)) = image2[k][i][j];
						if(actWindowAndGap2[k](i-(centre_pix_x-(int)(window_size/2)-max_gap),j-(centre_pix_y-(int)(window_size/2)-max_gap)) == 0)
						{
							output[0][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
							output[1][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
							output[2][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = 0;
							break_=1;
							break;
						}
					}
/* if calculation aborted: stop actual window calc.      */
			if (break_==1)
			{
				delete[] actWindowAndGap1;
				delete[] actWindowAndGap2;
				continue;
			}

/*********************************************************\
 * calculate cross-correlation for each possible gap     *
 * find gap with best cross-correlation                  *
\*********************************************************/
//--------------------------------------------------------
// define variables for identifying best cross-correlation
//--------------------------------------------------------
/* actualizing calculated pixels                         */
			long var3=*counter_used;
			*counter_used=var3+1;
/* array for best gap                                    */
			int best_local_gap[2];
			best_local_gap[0] = 0;
			best_local_gap[1] = 0;
/* max coefficient of determination of cross-correlation */
			double r_best = 0;
/* mean of all moving window values of image 1           */
			long double mean_x=0;
/* moving window without gap                             */
			MatrixXd* actWindow1 = new MatrixXd[bandsDataset1];
//--------------------------------------------------------
// calculate mean of all pixels values of moving window
//--------------------------------------------------------
			for (int k=0; k<bandsDataset1; k++)
			{
/* fill actual moving window                             */
				actWindow1[k] = actWindowAndGap1[k].block(max_gap,max_gap,window_size,window_size);
/* sum of moving window means of all bands               */
				mean_x += actWindow1[k].mean();
			}
/* mean of all band sums                                 */
			mean_x /= bandsDataset1;
//--------------------------------------------------------
// check windows of all possible gaps of image 2
//--------------------------------------------------------
			for (int gap_x=-max_gap; gap_x<=max_gap; gap_x++)
			{
				for (int gap_y=-max_gap; gap_y<=max_gap; gap_y++)
				{
/* define variables for statistic                        */
					long double mean_y=0, s_xy=0, s_xx=0, s_yy=0;
/* moving window of image 2 shift by actual gap          */
					MatrixXd* actWindow2 = new MatrixXd[bandsDataset2];
//--------------------------------------------------------
// calculate mean of all pixels values of moving window
//--------------------------------------------------------
					for (int k=0; k<bandsDataset2; k++)
					{
/* fill actual moving window of image 2                  */
						actWindow2[k] = actWindowAndGap2[k].block(max_gap+gap_x,max_gap+gap_y,window_size,window_size);
/* sum of moving window means of all bands               */
						mean_y += actWindow2[k].mean();
					}
/* mean of all band sums                                 */
					mean_y /= bandsDataset2;
//--------------------------------------------------------
// calculate cross-correlation of actual combination of shifted moving windows
//--------------------------------------------------------
					for (int window_z=0; window_z<bandsDataset2; window_z++)
					{
						for (int window_x=0; window_x<window_size; window_x++)
						{
							for (int window_y=0; window_y<window_size; window_y++)
							{
/* variances and covariances (temporal)                  */
								s_xx += ((actWindow1[window_z](window_x,window_y)-mean_x) * (actWindow1[window_z](window_x,window_y)-mean_x));
								s_xy += ((actWindow1[window_z](window_x,window_y)-mean_x) * (actWindow2[window_z](window_x,window_y)-mean_y));
								s_yy += ((actWindow2[window_z](window_x,window_y)-mean_y) * (actWindow2[window_z](window_x,window_y)-mean_y));
							}
						}
					}
/* variances and covariances                             */
					s_xx /= (window_size*window_size*bandsDataset1-1);
					s_xy /= (window_size*window_size*bandsDataset1-1);
					s_yy /= (window_size*window_size*bandsDataset1-1);
/* cross-correlation                                     */
					long double r_xy = s_xy/(sqrt(s_xx)*sqrt(s_yy));
					//std::cout << gap_x << " " << gap_y << " " << " " << s_xx << " " << s_yy << " " << s_xy << " " << r_xy << std::endl;
//--------------------------------------------------------
// check, if cross-correlation (coefficient of determination) is the highest
//--------------------------------------------------------
					if (r_best<=r_xy)
					{
/* if yes, remember as highest cross-correlation         */
						best_local_gap[0] = gap_x;
						best_local_gap[1] = gap_y;
						r_best = r_xy;
					}
/* clear memory                                          */
					delete[] actWindow2;
				}
			}
//--------------------------------------------------------
// extract gap with highest coefficient of determination
//--------------------------------------------------------
			output[0][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = best_local_gap[0];
			output[1][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = best_local_gap[1];
			output[2][centre_pix_x-((int)((window_size+1)/2)+max_gap-1)][centre_pix_y-((int)((window_size+1)/2)+max_gap-1)] = r_best;
/* clear memory                                          */
			delete[] actWindow1;
/* clear memory                                          */
			delete[] actWindowAndGap1;
			delete[] actWindowAndGap2;
/* user friendly text output                             */
		}
/*********************************************************\
 * END OF FUNCTION                                       *
\*********************************************************/
	}
/*********************************************************\
 * CONVERSION CONSTRUCTOR                                *
 * Input:                                                *
 * - image1: array of 1st multiband intersection image   *
 * - image2: array of 2nd multiband intersection image   *
 * - counter_abs: counter for theor. calculated pixels   *
 * - counter_used: counter for calculated pixels         *
 * - poDataset1: Image 1 (Metadata)                      *
 * - poDataset2: Image 2 (Metadata)                      *
 * - no_pix_x: number overlapping pixels in x-direction  *
 * - no_pix_y: number overlapping pixels in y-direction  *
 * - window_size: size of sliding window                 *
 * - max_gap: maximal tested gap of sliding window       *
 * - output: array for mean gap in x- and y-direction    *
 * - output2: array for gap of each intersection pixel   *
\*********************************************************/
	Parallelizing(float ***image1,float ***image2,float ***image_invar,int value_invar,long *counter_abs,long *counter_used,int bandsDataset1,int bandsDataset2,int no_pix_x,int no_pix_y,int window_size,int max_gap, float ***output)
		: image1(image1),image2(image2),image_invar(image_invar),value_invar(value_invar),counter_abs(counter_abs),counter_used(counter_used),bandsDataset1(bandsDataset1),bandsDataset2(bandsDataset2),no_pix_x(no_pix_x),no_pix_y(no_pix_y),window_size(window_size),max_gap(max_gap),output(output){}
};

GDALDataset* openImagefile(char *fname)
/*********************************************************\
 * FUNCTION 'openImagefile'                              *
 * Goal:                                                 *
 * - open an image dataset (read only)                   *
 * Input:                                                *
 * - path and name of image dataset                      *
 * Output:                                               *
 * - image dataset                                       *
\*********************************************************/
{
	GDALDataset *dataset;
	dataset = (GDALDataset *)GDALOpen(fname, GA_ReadOnly);
	if( dataset == NULL )
		fprintf(stderr,"Error: Image file not found.\n");
	return dataset;
}

OGRGeometry* buildGeometry (double ulx, double uly, double lrx, double lry, OGRSpatialReference * ref)
/*********************************************************\
 * FUNCTION 'buildGeometry'                              *
 * Goal:                                                 *
 * - create vector geometry based on given vertices      *
 * Input:                                                *
 * - given vertices                                      *
 * - coordinate system (call-by-reference)               *
 * Output:                                               *
 * - coordinate system (call-by-reference)               *
 * - vector geometry                                     *
\*********************************************************/
{
	char ulx_char[24], uly_char[24], lrx_char[24], lry_char[24];
/* changing variable format of ul and lr coordinates     */
	gcvt(ulx, 15, ulx_char);
	gcvt(uly, 15, uly_char);
	gcvt(lrx, 15, lrx_char);
	gcvt(lry, 15, lry_char);
/* changing variable format of ul and lr coordinates     */
	std::stringstream Wkt_temp;
	Wkt_temp << "POLYGON((" <<
		ulx_char << " " << uly_char << ", " << 
		lrx_char << " " << uly_char << ", " << 
		lrx_char << " " << lry_char << ", " <<
		ulx_char << " " << lry_char << ", " <<
		ulx_char << " " << uly_char << "))" << std::endl;
/* changing variable format of ul and lr coordinates     */
	std::string Wkt_temp_b = Wkt_temp.str();
/* changing variable format of ul and lr coordinates     */
	const char* Wkt_temp_c = Wkt_temp_b.c_str();
/* changing variable format of ul and lr coordinates     */
	char* Wkt = (char*) Wkt_temp_c;
/* creating a blank intersection image with calc. coord. */
	OGRGeometry * geom = OGRGeometryFactory::createGeometry(wkbPolygon);
	OGRErr ogrerr = OGRGeometryFactory::createFromWkt(&Wkt, ref, &geom);
	return geom;
}

OGREnvelope* getOverlap(double ulx1, double uly1, double lrx1, double lry1, double ulx2, double uly2, double lrx2, double lry2)
/*********************************************************\
 * FUNCTION 'getOverlap'                                 *
 * Goal:                                                 *
 * - create vector geometry of overlap area              *
 * Input:                                                *
 * - given vertices of two partially overlapping images  *
 * Output:                                               *
 * - vector geometry                                     *
\*********************************************************/
{
	OGREnvelope * poE = new OGREnvelope();
	poE->MinX = (ulx1 >= ulx2) ? ulx1 : ulx2;
	poE->MaxY = (uly1 <= uly2) ? uly1 : uly2;
	poE->MaxX = (lrx1 <= lrx2) ? lrx1 : lrx2;
	poE->MinY = (lry1 >= lry2) ? lry1 : lry2;
	return poE;
}

float*** arrayOfImage3D(GDALDataset * poDataset, OGREnvelope * poE, double refGeocode[6])
/*********************************************************\
 * FUNCTION 'arrayOfImage3D'                             *
 * Goal:                                                 *
 * - convert multi band image into 3D array              *
 * Input:                                                *
 * - given vertices of two partially overlapping images  *
 * Output:                                               *
 * - vector geometry                                     *
\*********************************************************/
{
	float ***image = 0;
	image = new float **[poDataset->GetRasterCount()];
	for (int i=0; i<poDataset->GetRasterCount(); i++)
	{
		image[i] = new float *[(int)((poE->MaxX-poE->MinX)/refGeocode[1])];
		for (int j=0; j<(int)((poE->MaxX-poE->MinX)/refGeocode[1]); j++)
			image[i][j] = new float [(int)((poE->MaxY-poE->MinY)/(-refGeocode[5]))];
	}
	//for (int i=0; i<poDataset->GetRasterCount(); i++)
	//{
	//	image[i] = new float *[(int)((poE->MaxY-poE->MinY)/-refGeocode[5])];
	//	for (int j=0; j<(int)((poE->MaxY-poE->MinY)/(-refGeocode[5])); j++)
	//		image[i][j] = new float[(int)((poE->MaxX-poE->MinX)/refGeocode[1])];
	//}
//--------------------------------------------------------
//	defining temporal intersection matrix for one line of image
//--------------------------------------------------------
	for (int i=0; i<(int)((poE->MaxY-poE->MinY)/(-refGeocode[5])); i++) // for all y
	{
		float** inter_orig = new float*[poDataset->GetRasterCount()];
/* fill unique bands of one line with values of image 1  */
		for (int j=1; j<=poDataset->GetRasterCount(); j++)
		{
			GDALRasterBand * poBand = poDataset->GetRasterBand(j);
			inter_orig[j-1] = new float[(int)((poE->MaxX-poE->MinX)/refGeocode[1])];
			poBand->RasterIO(GF_Read,
				ld2l((poE->MinX-refGeocode[0])/refGeocode[1]),
				ld2l((refGeocode[3]-poE->MaxY)/-refGeocode[5])+i,
				(long)((poE->MaxX-poE->MinX)/refGeocode[1]),
				1,
				inter_orig[j-1],
				(long)((poE->MaxX-poE->MinX)/refGeocode[1]),
				1,
				GDT_Float32,
				0,
				0);
/* copy line in intersection matrix of image 1           */
			for (int k=0; k<(int)((poE->MaxX-poE->MinX)/refGeocode[1]); k++) // for all x
				//image[j-1][i][k] = inter_orig[j-1][k];
				image[j-1][k][i] = inter_orig[j-1][k];
/* clear memory                                          */
			delete[] inter_orig[j-1];
		}
/* clear memory                                          */
		delete[] inter_orig;
	}
	return image;
}

int main(int argc, char* argv[])
/*********************************************************\
 * FUNCTION: MAIN                                        *
 * Goal:                                                 *
 * - calculate cross correlation between two images      *
 *   - all or specific pixels                            *
 *   - for image overlap                                 *
 * Input:                                                *
 * - path of input images                                *
 * - size sliding window                                 *
 * - maximal controlled gap                              *
 * Output:                                               *
 * - pixel gap, depending on wished output format        *
 *   - mean gap in x- and y-direction                    *
 *   - image of gaps for each overlapping pixel          *
\*********************************************************/
{
/*********************************************************\
 * Register open source libraries GDAL and OGR           *
\*********************************************************/
	GDALAllRegister();
	OGRRegisterAll();

/*********************************************************\
 * Declare user variables                                *
 * - pathes and names of input and output images         *
 * - maximal gap and moving window size                  *
\*********************************************************/
	std::cout << std::endl << "Get user specifications";
	char *pszSrcFilename1, *pszSrcFilename2, *pszDstFilename, *pszInvarFilename;
	int max_gap, window_size, value_invar;

/*********************************************************\
 * Read text file with user definition                   *
\*********************************************************/
//--------------------------------------------------------
//  Check user input
//--------------------------------------------------------
	if (argc != 15)
	{
		fprintf(stderr,"\t...Error: wrong number of infomations in command line\n");
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
	for (int i=1; i<argc; i+=2)
	{
		char *value_char;
		if (i+1 != argc)
		{
			value_char = argv[i+1];
			if (strcmp(argv[i],"-in1") == 0)
/* Path and filename of image dataset 1                  */
			{
				pszSrcFilename1 = new char[strlen(value_char)+1];
				strcpy(pszSrcFilename1,value_char);
				std::cout << "\t...input image dataset 1" << std::endl;
			} else if (strcmp(argv[i],"-in2") == 0)
/* Path and filename of image dataset 2                  */
			{
				pszSrcFilename2 = new char[strlen(value_char)+1];
				strcpy(pszSrcFilename2,value_char);
				std::cout << "\t\t\t...input image dataset 2" << std::endl;
			} else if (strcmp(argv[i],"-in_invar") == 0)
/* Path and filename of pseudo invariant features        */
			{
				pszInvarFilename = new char[strlen(value_char)+1];
				strcpy(pszInvarFilename,value_char);
				std::cout << "\t\t\t...invariant pixels dataset" << std::endl;
			} else if (strcmp(argv[i],"-val_invar") == 0)
/* Value of pseudo invariant features                    */
			{
				char* end;
				value_invar = (int)strtol(value_char,&end,10);
				if(*end != '\0')
				{
					std::cerr << "\t\t\t...Error: value for invariant pixels wrong defined" << std::endl;
					std::cerr << "\t\t\t...ABORTION" << std::endl;
					return EXIT_FAILURE;
				}
				std::cout << "\t\t\t...value for invariant pixels" << std::endl;
			} else if (strcmp(argv[i],"-out") == 0)
/* Path and filename of output image dataset             */
			{
				pszDstFilename = new char[strlen(value_char)+1];
				strcpy(pszDstFilename,value_char);
				std::cout << "\t\t\t...output image dataset" << std::endl;
			} else if (strcmp(argv[i],"-window_size") == 0)
/* Moving window size                                    */
			{
				char* end;
				window_size = (int)strtol(value_char,&end,10);
				if(*end != '\0')
				{
					std::cerr << "\t\t\t...Error: moving window size wrong defined" << std::endl;
					std::cerr << "\t\t\t...ABORTION" << std::endl;
					return EXIT_FAILURE;
				}
				std::cout << "\t\t\t...moving window size" << std::endl << "\t\t\t   " << window_size << "x" << window_size << " pixels" << std::endl;
			} else if (strcmp(argv[i],"-max_gap") == 0)
/* Moving window size                                    */
			{
				char* end;
				max_gap = (int)strtol(value_char,&end,10);
				if(*end != '\0')
				{
					std::cerr << "\t\t\t...Error: maximal gap wrong defined" << std::endl;
					std::cerr << "\t\t\t...ABORTION" << std::endl;
					return EXIT_FAILURE;
				}
				std::cout << "\t\t\t...maximal gap" << std::endl << "\t\t\t   " << max_gap << " pixels" << std::endl;
			} else
			{
				std::cerr << "\t\t\t...Error: not enough or invalid arguments" << std::endl;
				std::cerr << "\t\t\t...ABORTION" << std::endl;
				return EXIT_FAILURE;
			}
		}
	}

/*********************************************************\
 * Load input image 1                                    *
\*********************************************************/
	std::cout << std::endl << "Image 1";
	GDALDataset *poSrcDataset1;
	poSrcDataset1 = openImagefile(pszSrcFilename1);
//--------------------------------------------------------
//  Checking, if dataset exist
//--------------------------------------------------------
	if( poSrcDataset1 == NULL )
	{
		std::cerr << "\t\t\t...Error: file not found" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
//  Checking information about projection
//--------------------------------------------------------
	if( poSrcDataset1->GetProjectionRef()  == NULL )
	{
		std::cerr << "\t\t\tError: no information about projection" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
//  Checking information about pixel size
//--------------------------------------------------------
	double adfGeoTransform1[6];
	if( poSrcDataset1->GetGeoTransform( adfGeoTransform1 ) != CE_None )
	{
		std::cerr << "\t\t\tError: no information about pixel size" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
// Getting and writing dataset information
//--------------------------------------------------------
	std::cout << "\t\t\t...loaded" << std::endl;
/* Type of image                                         */
	std::cout << "\t\t\t   Type: " << poSrcDataset1->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME) << std::endl;
/* no of pixels in row and column, no of bands           */
	std::cout << "\t\t\t   Size: " << poSrcDataset1->GetRasterXSize() << " x " << poSrcDataset1->GetRasterYSize() << " x " << poSrcDataset1->GetRasterCount() << std::endl;
/* pixel size                                            */
	std::cout << "\t\t\t   Pixel Size: " << adfGeoTransform1[1] << " x " << adfGeoTransform1[5] << std::endl;

/*********************************************************\
 * Load input image 2                                    *
\*********************************************************/
	std::cout << std::endl << "Image 2";
	GDALDataset *poSrcDataset2;
	poSrcDataset2 = openImagefile(pszSrcFilename2);
//--------------------------------------------------------
//  Checking, if dataset exist
//--------------------------------------------------------
	if( poSrcDataset2 == NULL )
	{
		std::cerr << "\t\t\t...Error: file not found" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
//  Checking information about projection
//--------------------------------------------------------
	if( poSrcDataset2->GetProjectionRef()  == NULL )
	{
		std::cerr << "\t\t\tError: no information about projection" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
//  Checking information about pixel size
//--------------------------------------------------------
	double adfGeoTransform2[6];
	if( poSrcDataset2->GetGeoTransform( adfGeoTransform2 ) != CE_None )
	{
		std::cerr << "\t\t\tError: no information about pixel size" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
// Getting and writing dataset information
//--------------------------------------------------------
	std::cout << "\t\t\t...loaded" << std::endl;
/* Type of image                                         */
	std::cout << "\t\t\t   Type: " << poSrcDataset2->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME) << std::endl;
/* no of pixels in row and column, no of bands           */
	std::cout << "\t\t\t   Size: " << poSrcDataset2->GetRasterXSize() << " x " << poSrcDataset2->GetRasterYSize() << " x " << poSrcDataset2->GetRasterCount() << std::endl;
/* pixel size                                            */
	std::cout << "\t\t\t   Pixel Size: " << adfGeoTransform2[1] << " x " << adfGeoTransform2[5] << std::endl;

/*********************************************************\
 * Load invariant pixels                                 *
\*********************************************************/
	std::cout << std::endl << "Invariant pixels";
	GDALDataset *poInvarDataset;
	poInvarDataset = openImagefile(pszInvarFilename);
//--------------------------------------------------------
//  Checking, if dataset exist
//--------------------------------------------------------
	if( poInvarDataset == NULL )
	{
		std::cerr << "\t...Error: file not found" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
//  Checking information about projection
//--------------------------------------------------------
	if( poInvarDataset->GetProjectionRef()  == NULL )
	{
		std::cerr << "\tError: no information about projection" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
//  Checking information about pixel size
//--------------------------------------------------------
	double adfGeoTransform_invar[6];
	if( poInvarDataset->GetGeoTransform( adfGeoTransform_invar ) != CE_None )
	{
		std::cerr << "\tError: no information about pixel size" << std::endl;
		std::cerr << "\t\t\t...ABORTION" << std::endl;
		return EXIT_FAILURE;
	}
//--------------------------------------------------------
// Getting and writing dataset information
//--------------------------------------------------------
	std::cout << "\t...loaded" << std::endl;
/* Type of image                                         */
	std::cout << "\t\t\t   Type: " << poInvarDataset->GetDriver()->GetMetadataItem(GDAL_DMD_LONGNAME) << std::endl;
/* no of pixels in row and column, no of bands           */
	std::cout << "\t\t\t   Size: " << poInvarDataset->GetRasterXSize() << " x " << poInvarDataset->GetRasterYSize() << " x " << poInvarDataset->GetRasterCount() << std::endl;
/* pixel size                                            */
	std::cout << "\t\t\t   Pixel Size: " << adfGeoTransform_invar[1] << " x " << adfGeoTransform_invar[5] << std::endl;

/*********************************************************\
 * Extract images positions                              *
\*********************************************************/
	long double image1_ulx = adfGeoTransform1[0];
	long double image1_uly = adfGeoTransform1[3];
	long double image1_lrx = adfGeoTransform1[0] + adfGeoTransform1[1]*poSrcDataset1->GetRasterXSize();
	long double image1_lry = adfGeoTransform1[3] + adfGeoTransform1[5]*poSrcDataset1->GetRasterYSize();
	long double image2_ulx = adfGeoTransform2[0];
	long double image2_uly = adfGeoTransform2[3];
	long double image2_lrx = adfGeoTransform2[0] + adfGeoTransform2[1]*poSrcDataset2->GetRasterXSize();
	long double image2_lry = adfGeoTransform2[3] + adfGeoTransform2[5]*poSrcDataset2->GetRasterYSize();
	long double imageInvar_ulx = adfGeoTransform_invar[0];
	long double imageInvar_uly = adfGeoTransform_invar[3];
	long double imageInvar_lrx = adfGeoTransform_invar[0] + adfGeoTransform_invar[1]*poInvarDataset->GetRasterXSize();
	long double imageInvar_lry = adfGeoTransform_invar[3] + adfGeoTransform_invar[5]*poInvarDataset->GetRasterYSize();

/*********************************************************\
 * Create intersection arrays                            *
\*********************************************************/
	std::cout << std::endl << "Intersection matrizes";
	OGRSpatialReference * ref1 = new OGRSpatialReference(NULL);
	OGRGeometry * geom1 = buildGeometry(image1_ulx,image1_uly,image1_lrx,image1_lry,ref1);
	OGRSpatialReference * ref2 = new OGRSpatialReference(NULL);
	OGRGeometry * geom2 = buildGeometry(image2_ulx,image2_uly,image2_lrx,image2_lry,ref2);

/*********************************************************\
 * Find overlap area of neighboured images               *
\*********************************************************/
	OGRBoolean is_inter = geom1->Intersects(geom2);
	if (!is_inter)
	{
		std::cerr << "\tError: no overlap area." << std::endl;
		return EXIT_FAILURE;
	}
	OGREnvelope * poE = getOverlap(image1_ulx,image1_uly,image1_lrx,image1_lry,image2_ulx,image2_uly,image2_lrx,image2_lry);
	std::cout << "\t...defined" << std::endl;

/*********************************************************\
 * Define intersection matrices                          *
 * IMPORTANT: SAME GEOMETRIC RESOLUTION                  *
\*********************************************************/
/* image 1                                               */
	float *** image1 = arrayOfImage3D(poSrcDataset1, poE, adfGeoTransform1);
	int bandsDataset1 = poSrcDataset1->GetRasterCount();
	std::cout << "\t\t\t...image 1 read" << std::endl;
/* image 2                                               */
	float *** image2 = arrayOfImage3D(poSrcDataset2, poE, adfGeoTransform2);
	int bandsDataset2 = poSrcDataset2->GetRasterCount();
	std::cout << "\t\t\t...image 2 read" << std::endl;
/* invar pixels image                                    */
	float *** image_invar = arrayOfImage3D(poInvarDataset, poE, adfGeoTransform_invar);
	int bandsDataset_invar = poInvarDataset->GetRasterCount();
	std::cout << "\t\t\t...invar pixels read" << std::endl;

/*********************************************************\
 * Calculate size of intersecton subset                  *
 * - Use subset of image 1                               *
\*********************************************************/
	int no_pix_x = (int)((poE->MaxX-poE->MinX)/  adfGeoTransform1[1]) -(int)((window_size+1)/2)-max_gap-((int)((window_size+1)/2)+max_gap-1)+1;
	int no_pix_y = (int)((poE->MaxY-poE->MinY)/(-adfGeoTransform1[5]))-(int)((window_size+1)/2)-max_gap-((int)((window_size+1)/2)+max_gap-1)+1;

/*********************************************************\
 * Define output image                                   *
 * - 2 bands                                             *
 * - x- and y-gap for each overlapping pixel             *
\*********************************************************/
	std::cout << std::endl << "Output image";
/* image                                                 */
	GDALDriver *poDriver;
	const char *pszFormat = "GTiff";
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	GDALDataset *poDstDS;
/* characteristica                                       */
	poDstDS = poDriver->Create(pszDstFilename,no_pix_x,no_pix_y,3,GDT_Float32,NULL);
	double adfDstGeoTransform[6] =
	{
		poE->MinX+((int)((window_size+1)/2)-1+max_gap)*adfGeoTransform1[1],
		adfGeoTransform1[1],
		0,
		poE->MaxY+((int)((window_size+1)/2)-1+max_gap)*adfGeoTransform1[5],
		0,
		adfGeoTransform1[5]
	};
	poDstDS->SetGeoTransform(adfDstGeoTransform);
	OGRSpatialReference *refDstDS = ref1->Clone();
	poDstDS->SetProjection(poSrcDataset1->GetProjectionRef());
	poDstDS->SetMetadataItem( "SourceAgency", "Organizacion de las Naciones Unidas para la Alimentacion y la Agricultura");
	std::cout << "\t\t...Metadata defined" << std::endl;

/*********************************************************\
 * Define output array                                   *
 * - 3 bands                                             *
 * - x- and y-gap for each overlapping pixel             *
 * - coefficient of determination                        *
\*********************************************************/
	float ***output;
	output = new float **[3];
	for(int i=0; i<3; i++)
	{
		output[i] = new float *[no_pix_x];
		for (int j=0; j<no_pix_x; j++)
			output[i][j] = new float [no_pix_y];
	}

/*********************************************************\
 * Parallel Processing                                   *
 * - pixelwise calculation                               *
\*********************************************************/
	long counter_abs=0, counter_used=0;
	parallel_for(blocked_range<size_t>(0,no_pix_x*no_pix_y),Parallelizing(image1,image2,image_invar,value_invar,&counter_abs,&counter_used,bandsDataset1,bandsDataset2,no_pix_x,no_pix_y,window_size,max_gap,output));

/*********************************************************\
 * Save output image                                     *
\*********************************************************/
	std::cout << std::endl << "Save output image\t";
    int   nXSize = poDstDS->GetRasterXSize();
	int   nYSize = poDstDS->GetRasterYSize();
	float *bandArray;
	for (int i=0; i<3; i++)
	{
		int y_start = 0;
		while(y_start<no_pix_y)
		{
			int x_start = 0;
			int x_size = no_pix_x;
			int y_size = (y_start + (int)(1000000/no_pix_x) < no_pix_y) ? (int)(1000000/no_pix_x) : no_pix_y-y_start;
			float *bandArray;
			bandArray = (float *) CPLMalloc(sizeof(float)*(long)(x_size*y_size));
			for (int j=0; j<x_size; j++)
			{
				for (int k=0; k<y_size; k++)
				{
					bandArray[(long)(k*x_size+j)] = (float)output[i][x_start+j][y_start+k];
				}
			}
			poDstDS->GetRasterBand(i+1)->RasterIO(GF_Write,x_start,y_start,x_size,y_size,bandArray,x_size,y_size,GDT_Float32,0,0);
			y_start += (int)(1000000/no_pix_x);
			CPLFree(bandArray);
		}
	}
/* user friendly text output                             */
	std::cout << "...finished" << std::endl;

/*********************************************************\
 * clear memory                                          *
\*********************************************************/
/* output                                                */
	for(int i=0; i<3; i++)
	{
		for (int j=0; j<no_pix_x; j++)
			delete[] output[i][j];
		delete[] output[i];
	}
	delete[] output;

/* image 1                                               */
	for (int i=0; i<poSrcDataset1->GetRasterCount(); i++)
	{
		for (int j=0; j<(int)((poE->MaxX-poE->MinX)/adfGeoTransform1[1]); j++)
			delete[] image1[i][j];
		delete[] image1[i];
	}
	delete[] image1;
/* image 2                                               */
	for (int i=0; i<poSrcDataset1->GetRasterCount(); i++)
	{
		for (int j=0; j<(int)((poE->MaxX-poE->MinX)/adfGeoTransform2[1]); j++)
			delete[] image2[i][j];
		delete[] image2[i];
	}
	delete[] image2;
/* invar pixels image                                    */
	for (int i=0; i<poInvarDataset->GetRasterCount(); i++)
	{
		for (int j=0; j<(int)((poE->MaxX-poE->MinX)/adfGeoTransform_invar[1]); j++)
			delete[] image_invar[i][j];
		delete[] image_invar[i];
	}
	delete[] image_invar;

/*********************************************************\
 * closing datasets                                      *
\*********************************************************/
	GDALClose((GDALDatasetH) poSrcDataset1);
	GDALClose((GDALDatasetH) poSrcDataset2);
	GDALClose((GDALDatasetH) poDstDS);

/*********************************************************\
 * end of program                                        *
\*********************************************************/
	return 0;
}

long ld2l(long double f)
/*********************************************************\
 * FUNCTION: ld2l                                        *
 * Goal:                                                 *
 * - round single long double value                      *
 * Input:                                                *
 * - long double value                                   *
 * Output:                                               *
 * - long value                                          *
\*********************************************************/
{return f<0 ? f-.5 : f+.5;}
