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
#ifndef itkWaveletFrequencyInverse_hxx
#define itkWaveletFrequencyInverse_hxx
#include <itkWaveletFrequencyInverse.h>
#include <itkCastImageFilter.h>
#include <itkImage.h>
#include <algorithm>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkFrequencyExpandImageFilter.h>
#include <itkFrequencyExpandViaInverseFFTImageFilter.h>
#include <itkChangeInformationImageFilter.h>
namespace itk
{
template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
WaveletFrequencyInverse< TInputImage, TOutputImage, TWaveletFilterBank>
::WaveletFrequencyInverse()
  : m_Levels(1),
  m_HighPassSubBands(1),
  m_ScaleFactor(2)
{
  this->SetNumberOfRequiredOutputs(1);
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
std::pair<unsigned int, unsigned int>
WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::InputIndexToLevelBand(unsigned int linear_index)
{
  if (linear_index > this->m_TotalInputs - 1 || linear_index < 0)
    itkExceptionMacro(<< "Failed converting liner index " << linear_index
                      << " to Level,Band pair : out of bounds");
  // Low pass (band = 0).
  if (linear_index == 0 )
    return std::make_pair(this->m_Levels, 0);

  unsigned int band = (linear_index - 1) % this->m_HighPassSubBands;
  band = band + 1;
  // note integer division ahead.
  unsigned int level = (linear_index - band) / this->m_HighPassSubBands + 1;
  return std::make_pair(level, band);
};

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetLevels(unsigned int n)
{
  unsigned int current_inputs = 1 + this->m_Levels * this->m_HighPassSubBands;

  if ( this->m_TotalInputs == current_inputs && this->m_Levels == n )
    {
    return;
    }

  this->m_Levels = n;
  this->m_TotalInputs = 1 + n * this->m_HighPassSubBands;

  this->SetNumberOfRequiredInputs( this->m_TotalInputs );
  this->Modified();

  this->SetNthOutput(0, this->MakeOutput(0));
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetHighPassSubBands(unsigned int k)
{
  if ( this->m_HighPassSubBands == k )
    {
    return;
    }
  this->m_HighPassSubBands = k;
  // Trigger setting new number of inputs avoiding code duplication
  this->SetLevels(this->m_Levels);
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void
WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetInputs(const std::vector<InputImagePointer> & inputs)
{
  if (inputs.size() != this->m_TotalInputs)
    itkExceptionMacro(<< "Error seting inputs in inverse wavelet. Wrong vector size: "
                      << inputs.size() << " .According to number of levels and bands it should be: " << m_TotalInputs);

  for( unsigned int nin = 0; nin < this->m_TotalInputs; ++nin)
    {
    if( this->GetInput(nin) != inputs[nin])
      this->SetNthInput(nin, inputs[nin]);
    }
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void
WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetInputLowPass(const InputImagePointer & input_low_pass)
{
  this->SetNthInput(0, input_low_pass);
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void
WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::SetInputsHighPass(const std::vector<InputImagePointer> & inputs)
{
  if (inputs.size() != this->m_TotalInputs - 1)
    itkExceptionMacro(<< "Error seting inputs in inverse wavelet. Wrong vector size: "
                      << inputs.size() << " .According to number of levels and bands it should be: " << m_TotalInputs
      - 1);

  for( unsigned int nin = 0; nin < this->m_TotalInputs - 1; ++nin)
    {
    this->SetNthInput(nin + 1, inputs[nin]);
    }
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank >
void WaveletFrequencyInverse< TInputImage, TOutputImage, TWaveletFilterBank>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent
     << " Levels: " << this->m_Levels
     << " HighPassSubBands: " << this->m_HighPassSubBands
     << " TotalInputs: " << this->m_TotalInputs
     << std::endl;
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  // Check  all inputs exist.
  for (unsigned int n_input = 0; n_input < this->m_TotalInputs; ++n_input)
    {
    if (!this->GetInput(n_input))
      itkExceptionMacro(<< "Input: " << n_input << " has not been set");
    }

  // We know inputIndex = 1 has the same size than output. Use it.
  InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput(1));
  const typename InputImageType::PointType & inputOrigin =
    inputPtr->GetOrigin();
  const typename InputImageType::SpacingType & inputSpacing =
    inputPtr->GetSpacing();
  const typename InputImageType::DirectionType & inputDirection =
    inputPtr->GetDirection();
  const typename InputImageType::SizeType & inputSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  const typename InputImageType::IndexType & inputStartIndex =
    inputPtr->GetLargestPossibleRegion().GetIndex();

  OutputImagePointer outputPtr;
  outputPtr = this->GetOutput(0);

  typename OutputImageType::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize(inputSize);
  outputLargestPossibleRegion.SetIndex(inputStartIndex);

  outputPtr->SetLargestPossibleRegion(outputLargestPossibleRegion);
  outputPtr->SetOrigin(inputOrigin);
  outputPtr->SetSpacing(inputSpacing);
  outputPtr->SetDirection(inputDirection);
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateOutputRequestedRegion(DataObject *refOutput)
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputRequestedRegion(refOutput);
}

template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse<TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // compute baseIndex and baseSize
  typedef typename OutputImageType::SizeType   SizeType;
  typedef typename OutputImageType::IndexType  IndexType;
  typedef typename OutputImageType::RegionType RegionType;

  OutputImagePointer outputPtr = this->GetOutput(0);
  std::cout << "Output 0 :" <<  outputPtr->GetRequestedRegion() << '\n';

  IndexType  inputIndex;
  SizeType   inputSize;
  RegionType inputRegion;
  SizeType   baseSize  = outputPtr->GetRequestedRegion().GetSize();
  IndexType  baseIndex = outputPtr->GetRequestedRegion().GetIndex();
  RegionType baseRegion;
  baseRegion.SetIndex(baseIndex);
  baseRegion.SetSize(baseSize);
  inputRegion = baseRegion;

  for (unsigned int level = 0; level < this->m_Levels; ++level)
    {
    for (unsigned int band = 0; band < this->m_HighPassSubBands; ++band)
      {
      unsigned int n_input = 1 + level * this->m_HighPassSubBands + band;
      if (!this->GetInput(n_input))
        continue;
      InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput(n_input));
      // make sure the region is within the largest possible region
      std::cout << "GenerateINPUT : " << n_input <<  inputPtr->GetLargestPossibleRegion() << '\n';
      inputRegion.Crop(inputPtr->GetLargestPossibleRegion());
      // set the requested region
      inputPtr->SetRequestedRegion(inputRegion);
      }

    /******* Update base region for next level *********/
    unsigned int scaleFactorPerLevel = std::pow(static_cast<double>(this->m_ScaleFactor), static_cast<int>(level + 1));
    for (unsigned int idim = 0; idim < TInputImage::ImageDimension; idim++)
      {
      // inputIndex[idim] = baseIndex[idim] * scaleFactorPerLevel;
      // inputSize[idim] = baseSize[idim] * scaleFactorPerLevel;
      // Index by half.
      inputIndex[idim] = static_cast<IndexValueType>(
          std::ceil(static_cast<double>(baseIndex[idim]) / scaleFactorPerLevel));
      // Size by half
      inputSize[idim] = static_cast<SizeValueType>(
          std::floor(static_cast<double>(baseSize[idim]) / scaleFactorPerLevel));
      if (inputSize[idim] < 1)
        itkExceptionMacro(
            << "Failure at level: " << level
            << " in forward wavelet, going to negative image size. Too many levels for input image size.");
      }

    // Update Base Region for next levels.
    inputRegion.SetIndex(inputIndex);
    inputRegion.SetSize(inputSize);

    // Set low pass input
    if(level == this->m_Levels - 1)
      {
      unsigned int n_input = 0;
      if (!this->GetInput(n_input))
        continue;
      InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput(n_input));
      std::cout << "GenerateINPUT.LowPass Input : " << n_input <<  inputPtr->GetLargestPossibleRegion() << '\n';
      inputRegion.Crop(inputPtr->GetLargestPossibleRegion());
      inputPtr->SetRequestedRegion(inputRegion);
      }
    }
}

// ITK forward implementation: Freq Domain
//    - HPs (lv1 wavelet coef)
// I -             - HPs (lv2 wavelet coef)
//    - LP * Down -
//                 - LP * Down
// Where Down is a downsample.
template<typename TInputImage, typename TOutputImage, typename TWaveletFilterBank>
void WaveletFrequencyInverse< TInputImage, TOutputImage, TWaveletFilterBank>
::GenerateData()
{
  this->AllocateOutputs();
  InputImageConstPointer low_pass = this->GetInput(0);

  typedef itk::CastImageFilter<InputImageType, OutputImageType> CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(low_pass);
  castFilter->Update();
  OutputImagePointer low_pass_per_level = castFilter->GetOutput();
  // DEBUG info
    {
    typename OutputImageType::IndexType zeroIndex;
    zeroIndex.Fill(0);
    std::cout << "low_pass_per_level: input:" << low_pass_per_level->GetLargestPossibleRegion() << "DC component:"
              << low_pass_per_level->GetPixel(zeroIndex) << '\n';
    }

  for (int level = this->m_Levels - 1; level > -1; --level)
    {
    std::cout << "LEVEL: " << level << std::endl;
    typedef itk::MultiplyImageFilter<OutputImageType> MultiplyFilterType;
    /******** Upsample LowPass ********/
// #ifdef ITK_VISUALIZE_TESTS
//   typedef itk::Image<float, ImageDimension> SpatialImageType;
//   typedef itk::InverseFFTImageFilter<OutputImageType, SpatialImageType> InverseFFTFilterType;
//   typename InverseFFTFilterType::Pointer inverseFFT = InverseFFTFilterType::New();
//   inverseFFT->SetInput(low_pass_per_level);
//   inverseFFT->Update();
//   Testing::ViewImage(inverseFFT->GetOutput(), "Approx before upsampling");
// #endif
    // typedef itk::FrequencyExpandViaInverseFFTImageFilter<OutputImageType> ExpandFilterType;
    typedef itk::FrequencyExpandImageFilter<OutputImageType> ExpandFilterType;
    typename ExpandFilterType::Pointer upsampleFilter = ExpandFilterType::New();
    upsampleFilter->SetInput(low_pass_per_level);
    upsampleFilter->SetExpandFactors(this->m_ScaleFactor);
    upsampleFilter->Update();
// #ifdef ITK_VISUALIZE_TESTS
//     inverseFFT->SetInput(upsampleFilter->GetOutput());
//     inverseFFT->Update();
//     Testing::ViewImage(inverseFFT->GetOutput(), "Approx after upsampling");
// #endif

    typename MultiplyFilterType::Pointer multiplyUpsampleCorrection = MultiplyFilterType::New();
    multiplyUpsampleCorrection->SetInput1(upsampleFilter->GetOutput());
    double expUpsampleCorrection = static_cast<double>(ImageDimension);
    multiplyUpsampleCorrection->SetConstant(std::pow(2.0, expUpsampleCorrection));
    multiplyUpsampleCorrection->InPlaceOn();
    multiplyUpsampleCorrection->Update();
// #ifdef ITK_VISUALIZE_TESTS
//     inverseFFT->SetInput(multiplyUpsampleCorrection->GetOutput());
//     inverseFFT->Update();
//     Testing::ViewImage(inverseFFT->GetOutput(), "Approx after upsampling correction");
// #endif
// Ignore modifications of origin and spacing of upsample filters.
    typedef itk::ChangeInformationImageFilter<OutputImageType> ChangeInformationFilterType;
    typename ChangeInformationFilterType::Pointer changeInfoFilter = ChangeInformationFilterType::New();
    // changeInfoFilter->SetInput(upsampleFilter->GetOutput());
    changeInfoFilter->SetInput(multiplyUpsampleCorrection->GetOutput());
    changeInfoFilter->ChangeDirectionOff();
    changeInfoFilter->ChangeRegionOff();
    changeInfoFilter->ChangeSpacingOn();
    changeInfoFilter->ChangeOriginOn();
    const typename OutputImageType::SpacingType new_spacing( 1.0 );
    changeInfoFilter->SetOutputSpacing(new_spacing);
    const typename OutputImageType::PointType new_origin( 0.0 );
    changeInfoFilter->SetOutputOrigin(new_origin);
    changeInfoFilter->Update();
    low_pass_per_level = changeInfoFilter->GetOutput();

    /******* Calculate FilterBank with the right size per level. *****/
      {
      // itkDebugMacro(<<"Low_pass_per_level: " << level << " Region:" << low_pass_per_level->GetLargestPossibleRegion() );
      typename OutputImageType::IndexType zeroIndex;
      zeroIndex.Fill(0);
      std::cout << "low_pass_per_level: " << level << " Region:" << low_pass_per_level->GetLargestPossibleRegion()
                << "DC component:" << low_pass_per_level->GetPixel(zeroIndex) << '\n';
      }
    typename WaveletFilterBankType::Pointer filterBank = WaveletFilterBankType::New();
    filterBank->SetHighPassSubBands(this->m_HighPassSubBands);
    filterBank->SetSize(low_pass_per_level->GetLargestPossibleRegion().GetSize() );
    filterBank->SetInverseBank(true);
    filterBank->Update();

    /******* LowPass band *****/
    typename MultiplyFilterType::Pointer multiplyLowPass = MultiplyFilterType::New();
    multiplyLowPass->SetInput1(filterBank->GetOutputLowPass());
    multiplyLowPass->SetInput2(low_pass_per_level);
    multiplyLowPass->InPlaceOn();
    multiplyLowPass->Update();
    low_pass_per_level = multiplyLowPass->GetOutput();
// #ifdef ITK_VISUALIZE_TESTS
//   inverseFFT->SetInput(low_pass_per_level);
//   inverseFFT->Update();
//   Testing::ViewImage(inverseFFT->GetOutput(), "Low Band x filter bank");
// #endif

    // if (level != static_cast<int>(this->m_Levels) - 1 )
    // {
    // #<{(|****** Dilation factor in the low_pass to have the same factor than input high bands of this level ****|)}>#
    // typename MultiplyFilterType::Pointer multiplyLowByReconstructLevelFactor = MultiplyFilterType::New();
    // multiplyLowByReconstructLevelFactor->SetInput1(multiplyLowPass->GetOutput());
    // double expLevelFactor = static_cast<double>(ImageDimension)/2.0;
    // multiplyLowByReconstructLevelFactor->SetConstant(std::pow(2.0, expLevelFactor));
    // multiplyLowByReconstructLevelFactor->InPlaceOn();
    // multiplyLowByReconstructLevelFactor->Update();
    // low_pass_per_level = multiplyLowByReconstructLevelFactor->GetOutput();
    // }

    /******* HighPass sub-bands *****/
    std::vector<OutputImagePointer> highPassMasks = filterBank->GetOutputsHighPassBands();
    // Store HighBands steps into high_pass_reconstruction image:
    OutputImagePointer reconstructed = OutputImageType::New();
    reconstructed->SetRegions(low_pass_per_level->GetLargestPossibleRegion());
    reconstructed->Allocate();
    reconstructed->FillBuffer(0);
    for(unsigned int band = 0; band < this->m_HighPassSubBands; ++band)
      {
      unsigned int n_input = 1 + level * this->m_HighPassSubBands + band;

      typename MultiplyFilterType::Pointer multiplyHighBandFilter = MultiplyFilterType::New();
      multiplyHighBandFilter->SetInput1(highPassMasks[band]);
      multiplyHighBandFilter->SetInput2(this->GetInput(n_input));
      multiplyHighBandFilter->InPlaceOn();
      multiplyHighBandFilter->Update();
// #ifdef ITK_VISUALIZE_TESTS
//   inverseFFT->SetInput(multiplyHighBandFilter->GetOutput());
//   inverseFFT->Update();
//   Testing::ViewImage(inverseFFT->GetOutput(), "High Pass x high wavelet filter");
// #endif

      /******* Band dilation factor for HighPass bands *****/
      //  2^(1/#bands) instead of Dyadic dilations.
      typename MultiplyFilterType::Pointer multiplyByReconstructionBandFactor = MultiplyFilterType::New();
      multiplyByReconstructionBandFactor->SetInput1(multiplyHighBandFilter->GetOutput());
      // double expFactorHigh = static_cast<int>(band + 1)/static_cast<double>(this->m_HighPassSubBands) * static_cast<double>(ImageDimension)/2.0;
      // double expFactorHigh = static_cast<double>(1)/static_cast<double>(this->m_HighPassSubBands) * static_cast<double>(ImageDimension)/2.0;
      // double expBandFactor = 0;// + static_cast<double>(ImageDimension)/2.0;
      // double expBandFactor = static_cast<int>(band)/static_cast<double>(this->m_HighPassSubBands) * static_cast<double>(ImageDimension)/2.0;
      double expBandFactor = 0;
      multiplyByReconstructionBandFactor->SetConstant(std::pow(2.0, expBandFactor));
      multiplyByReconstructionBandFactor->InPlaceOn();
      multiplyByReconstructionBandFactor->Update();

      typedef itk::AddImageFilter<OutputImageType> AddFilterType;
      typename AddFilterType::Pointer addFilter = AddFilterType::New();
      addFilter->SetInput1(reconstructed);
      addFilter->SetInput2(multiplyByReconstructionBandFactor->GetOutput());
      addFilter->InPlaceOn();
      addFilter->Update();
      reconstructed = addFilter->GetOutput();
// #ifdef ITK_VISUALIZE_TESTS
//   inverseFFT->SetInput(reconstructed);
//   inverseFFT->Update();
//   Testing::ViewImage(inverseFFT->GetOutput(), "reconstructed:Sum of High Bands ");
// #endif

      this->UpdateProgress(static_cast<float>(m_TotalInputs - n_input - 1)
                           / static_cast<float>(m_TotalInputs));
      }

    /******* Add low pass to the sum of high pass bands. *****/
    typedef itk::AddImageFilter<OutputImageType> AddFilterType;
    typename AddFilterType::Pointer addHighAndLow = AddFilterType::New();
    addHighAndLow->SetInput1(reconstructed.GetPointer()); // HighBands
    // addHighAndLow->SetInput2(multiplyLowByReconstructLevelFactor->GetOutput());
    addHighAndLow->SetInput2(low_pass_per_level);
    addHighAndLow->InPlaceOn();

    if (level == 0 /* Last level to compute */) // Graft Output
      {
      addHighAndLow->GraftOutput(this->GetOutput());
      addHighAndLow->Update();
// #ifdef ITK_VISUALIZE_TESTS
//   inverseFFT->SetInput(addHighAndLow->GetOutput());
//   inverseFFT->Update();
//   Testing::ViewImage(inverseFFT->GetOutput(), "add high and low without factor ");
// #endif
      this->GraftOutput(addHighAndLow->GetOutput());
      }
    else // Update low_pass
      {
      addHighAndLow->Update();
      /******* Dilation factor in the low_pass to have the same factor than input high bands of next level *****/
      typename MultiplyFilterType::Pointer multiplyLowByReconstructLevelFactor = MultiplyFilterType::New();
      multiplyLowByReconstructLevelFactor->SetInput1(addHighAndLow->GetOutput());
      // double expLevelFactor = static_cast<double>(ImageDimension)/2.0;
      double expLevelFactor = 0;
      multiplyLowByReconstructLevelFactor->SetConstant(std::pow(2.0, expLevelFactor));
      multiplyLowByReconstructLevelFactor->InPlaceOn();
      multiplyLowByReconstructLevelFactor->Update();
      low_pass_per_level = multiplyLowByReconstructLevelFactor->GetOutput();
      }
    }
}
} // end namespace itk
#endif
