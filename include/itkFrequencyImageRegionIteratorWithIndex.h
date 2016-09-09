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
#ifndef itkFrequencyImageRegionIteratorWithIndex_h
#define itkFrequencyImageRegionIteratorWithIndex_h

#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
/** \class FrequencyImageRegionIteratorWithIndex
 * \brief A multi-dimensional iterator templated over image type that walks
 * pixels within a region and is specialized to keep track of its image index
 * location.
 *
 * This class is a specialization of ImageRegionIteratorWithIndex,
 * adding method GetFrequencyIndex to give the frequency bins corresponding to image indices.
 * The frequency bins depends on the image size. The default assumes that the image to iterate over is
 * the output of a forward FFT filter, where the first index corresponds to 0 frequency, and Nyquist Frequencies are
 * in the middle, between positive and negative frequencies.
 *
 * This class can be specialized further to iterate over other frequency layouts, for example shifted images (where 0 frequency is in the middle of the image, and Nyquist are in the border).
 * The frequency layout is assumed to be: where fs is frequency sampling, or frequency spacing (1.0 by default).
 *
 * If N is even:
 * Nyquist frequency at index=N/2 is shared between + and - regions.
 * <------------positive f ---------------><------------negative f-------------->
 * 0(DC) fs/N 2*fs/N  ... (N/2 -1)*fs/N  fs/2   -(N/2-1)*fs/N  ...  -2*fs/N  -fs/N
 *
 * Example: Size 6:
 * +    |     -
 * ------------
 *      0       <-- DC Component (0 freq)
 * 1    |     5
 * 2    |     4
 *      3       <-- Shared between regions, unique Nyquist.
 *
 * If N is odd:
 * Nyquist frequency is not represented but there are simmetric largest frequencies at index=N/2, N/2 +1
 * <----------positive f ---------------><------------negative f----------------->
 * 0(DC) fs/N 2*fs/N  ...... fs/2*(N-1)/N    -fs/2*(N-1)/N  ...    -2*fs/N  -fs/N
 *
 * Example: Size 5:
 * +    |     -
 * ------------
 *      0       <-- DC Component (0 freq)
 * 1    |     4
 * 2    |     3 <-- Absolute Largest Frequency bins (+, -)
 *
 * Please see ImageRegionConstIteratorWithIndex for more information.
 * \sa ForwardFFTImageFilter
 *
 * \par MORE INFORMATION
 * For a complete description of the ITK Image Iterators and their API, please
 * see the Iterators chapter in the ITK Software Guide.  The ITK Software Guide
 * is available in print and as a free .pdf download from https://www.itk.org.
 *
 * \ingroup ImageIterators
 *
 * \sa ImageRegionIteratorWithIndex
 * \sa ImageConstIterator \sa ConditionalConstIterator
 * \sa ConstNeighborhoodIterator \sa ConstShapedNeighborhoodIterator
 * \sa ConstSliceIterator  \sa CorrespondenceDataStructureIterator
 * \sa FloodFilledFunctionConditionalConstIterator
 * \sa FloodFilledImageFunctionConditionalConstIterator
 * \sa FloodFilledImageFunctionConditionalIterator
 * \sa FloodFilledSpatialFunctionConditionalConstIterator
 * \sa FloodFilledSpatialFunctionConditionalIterator
 * \sa ImageConstIterator \sa ImageConstIteratorWithIndex
 * \sa ImageIterator \sa ImageIteratorWithIndex
 * \sa ImageLinearConstIteratorWithIndex  \sa ImageLinearIteratorWithIndex
 * \sa ImageRandomConstIteratorWithIndex  \sa ImageRandomIteratorWithIndex
 * \sa ImageRegionConstIterator \sa ImageRegionConstIteratorWithIndex
 * \sa ImageRegionExclusionConstIteratorWithIndex
 * \sa ImageRegionExclusionIteratorWithIndex
 * \sa ImageRegionIterator  \sa ImageRegionIteratorWithIndex
 * \sa ImageRegionReverseConstIterator  \sa ImageRegionReverseIterator
 * \sa ImageReverseConstIterator  \sa ImageReverseIterator
 * \sa ImageSliceConstIteratorWithIndex  \sa ImageSliceIteratorWithIndex
 * \sa NeighborhoodIterator \sa PathConstIterator  \sa PathIterator
 * \sa ShapedNeighborhoodIterator  \sa SliceIterator
 * \sa ImageConstIteratorWithIndex
 *
 * \ingroup IsotropicWavelets
 *
 */
template< typename TImage >
class FrequencyImageRegionIteratorWithIndex:public ImageRegionIteratorWithIndex< TImage >
{
public:
  /** Standard class typedefs. */
  typedef FrequencyImageRegionIteratorWithIndex  Self;
  typedef ImageRegionIteratorWithIndex< TImage > Superclass;

  /** Types inherited from the Superclass */
  typedef typename Superclass::IndexType             IndexType;
  typedef typename Superclass::SizeType              SizeType;
  typedef typename Superclass::OffsetType            OffsetType;
  typedef typename Superclass::RegionType            RegionType;
  typedef typename Superclass::ImageType             ImageType;
  typedef typename Superclass::PixelContainer        PixelContainer;
  typedef typename Superclass::PixelContainerPointer PixelContainerPointer;
  typedef typename Superclass::InternalPixelType     InternalPixelType;
  typedef typename Superclass::PixelType             PixelType;
  typedef typename Superclass::AccessorType          AccessorType;

  typedef typename ImageType::SpacingType            FrequencyIndexType;
  typedef typename ImageType::SpacingValueType       FrequencyIndexValueType;
  /** Default constructor. Needed since we provide a cast constructor. */
  FrequencyImageRegionIteratorWithIndex():
    ImageRegionIteratorWithIndex< TImage >()
  {
    this->InitIndices();
  };

  /** Constructor establishes an iterator to walk a particular image and a
   * particular region of that image. */
  FrequencyImageRegionIteratorWithIndex(TImage *ptr, const RegionType & region):
    ImageRegionIteratorWithIndex< TImage >(ptr, region)
  {
    this->InitIndices();
  };

  /** Constructor that can be used to cast from an ImageIterator to an
   * ImageRegionIteratorWithIndex. Many routines return an ImageIterator, but for a
   * particular task, you may want an ImageRegionIteratorWithIndex.  Rather than
   * provide overloaded APIs that return different types of Iterators, itk
   * returns ImageIterators and uses constructors to cast from an
   * ImageIterator to a ImageRegionIteratorWithIndex. */
  FrequencyImageRegionIteratorWithIndex(const ImageIteratorWithIndex< TImage > & it):
    ImageRegionIteratorWithIndex< TImage >(it)
  {
    this->InitIndices();
  };

  /** Note that this method is independent of the region in the constructor. **/
  FrequencyIndexType GetFrequencyIndex() const
  {
    FrequencyIndexType freq;
    freq.Fill(0);

    for (unsigned int dim = 0; dim < TImage::ImageDimension; dim++)
    {
      // Origin in a freq image should be always zero, we ignore it.
      // However spacing has a useful meaning: dF, delta(freq).
      // Freq: + (from 0 (DC) to largest f (Nyquist if even) in +dF intervals.
      if (this->m_PositionIndex[dim] <= m_HalfIndex[dim])
      {
        freq[dim] = this->m_Image->GetSpacing()[dim] *
          (this->m_PositionIndex[dim] - this->m_MinIndex[dim])
          / this->m_Image->GetLargestPossibleRegion().GetSize()[dim];
      }
      // Freq: -. From -Largest(Nyquist if even) to -dF frequency.
      else
      {
        freq[dim] = this->m_Image->GetSpacing()[dim] *
          (this->m_PositionIndex[dim] - (this->m_MaxIndex[dim] + 1))
          / this->m_Image->GetLargestPossibleRegion().GetSize()[dim];
      }
    }
    return freq;
  }

  FrequencyIndexValueType GetFrequencyModuloSquare() const
  {
    FrequencyIndexValueType w2(0);
    FrequencyIndexType w( this->GetFrequencyIndex() );

    for (unsigned int dim = 0; dim < TImage::ImageDimension; dim++)
    {
      w2 += w[dim] * w[dim];
    }
    return w2;
  }

  typename ImageType::IndexType GetHalfIndex() const
    {
    typename ImageType::IndexType half_index;
    for (unsigned int dim = 0; dim < ImageType::ImageDimension; dim++)
      {
      half_index[dim] = static_cast<typename ImageType::IndexValueType>(m_HalfIndex[dim]);
      }
    return half_index;
    };

protected:
  /** Calculate Nyquist frequency index (m_HalfIndex), and Min/Max indices from LargestPossibleRegion.
   * Called at constructors.  */
  void InitIndices()
    {
    for (unsigned int dim = 0; dim < ImageType::ImageDimension; dim++)
      {
      this->m_MinIndex[dim] =
        this->m_Image->GetLargestPossibleRegion().GetIndex()[dim];
      this->m_MaxIndex[dim] =
        this->m_Image->GetLargestPossibleRegion().GetIndex()[dim] +
        (this->m_Image->GetLargestPossibleRegion().GetSize()[dim] - 1 );
      this->m_HalfIndex[dim] = static_cast<FrequencyIndexValueType>(
        this->m_MinIndex[dim] +
        std::ceil( (this->m_MaxIndex[dim] - this->m_MinIndex[dim] ) / 2.0 )
        );
      }
    }
  /**
   * Index corresponding to the first highest frequency (Nyquist) after a FFT transform.
   * If the size of the image is even, the Nyquist frequency = fs/2 is unique and shared
   * between positive and negative frequencies.
   * If odd, Nyquist frequency is not represented, but there is still a largest frequency at this index
   * = fs/2 * (N-1)/N.
   */
  FrequencyIndexType m_HalfIndex;
  FrequencyIndexType m_MinIndex;
  FrequencyIndexType m_MaxIndex;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageRegionIteratorWithIndex.hxx"
#endif

#endif
