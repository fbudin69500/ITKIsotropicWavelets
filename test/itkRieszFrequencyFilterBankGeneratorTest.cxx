/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include <memory>
#include <string>
#include <cmath>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkForwardFFTImageFilter.h"
#include "itkRieszFrequencyFilterBankGenerator.h"
#include <itkComplexToRealImageFilter.h>
#include <itkImageRegionConstIterator.h>
//Visualize for dev/debug purposes. Set in cmake file. Require VTK
#if ITK_VISUALIZE_TESTS != 0
#include "itkView3DImage.h"
#endif

using namespace std;
using namespace itk;

int itkRieszFrequencyFilterBankGeneratorTest(int argc, char** argv)
{
    if( argc != 3 )
    {
        std::cerr << "Usage: " << argv[0] << " inputImage outputImage" << std::endl;
        return EXIT_FAILURE;
    }
    const string inputImage  = argv[1];
    const string outputImage = argv[2];

    const unsigned int dimension = 3;
    typedef double                           PixelType;
    typedef itk::Image<PixelType, dimension> ImageType;
    typedef itk::ImageFileReader<ImageType>  ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(inputImage);
    reader->Update();
    reader->UpdateLargestPossibleRegion();

    // Perform FFT on input image.
    typedef itk::ForwardFFTImageFilter<ImageType> FFTFilterType;
    typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
    fftFilter->SetInput(reader->GetOutput());
    fftFilter->Update();
    typedef typename FFTFilterType::OutputImageType ComplexImageType;

    // typedef itk::RieszFrequencyFunction<> FunctionType;
    typedef itk::RieszFrequencyFilterBankGenerator< ComplexImageType> RieszFilterBankType;
    typename RieszFilterBankType::Pointer filterBank = RieszFilterBankType::New();
    filterBank->SetSize(fftFilter->GetOutput()->GetLargestPossibleRegion().GetSize());
    filterBank->Update();

    //Get real part of complex image for visualization
    typedef itk::ComplexToRealImageFilter<ComplexImageType, ImageType> ComplexToRealFilter;
    typename ComplexToRealFilter::Pointer complexToRealFilter = ComplexToRealFilter::New();
    std::cout << "Real Part of ComplexImage:"<< std::endl;
    for (unsigned int dir = 0; dir < ImageType::ImageDimension; dir++)
    {
        std::cout << "Direction: " << dir + 1 << " / " << ImageType::ImageDimension << std::endl;
        complexToRealFilter->SetInput(filterBank->GetOutput(dir) );
        complexToRealFilter->Update();
#if ITK_VISUALIZE_TESTS != 0
        View3DImage(complexToRealFilter->GetOutput());
#endif
    }
    return EXIT_SUCCESS;
}
