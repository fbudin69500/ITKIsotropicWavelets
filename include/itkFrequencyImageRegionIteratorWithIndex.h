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
#include "itkFrequencyImageRegionConstIteratorWithIndex.h"

namespace itk
{
/** \class FrequencyImageRegionIteratorWithIndex

 * Iterator providing method GetFrequency() to retrieve the frequency associated to an index.
 * This value is related to the specific layout of frequencies from an image in the dual (frequency) space.
 * In this case, the layout corresponds to the output of a FastFourierTransform from FFTW library.
 * This class is a non-const version of FrequencyImageRegionConstIteratorWithIndex.
 *
 * \ingroup ImageIterators
 * \sa FrequencyImageRegionConstIteratorWithIndex
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
class FrequencyImageRegionIteratorWithIndex:
  public FrequencyImageRegionConstIteratorWithIndex<TImage>
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

  typedef typename ImageType::SpacingType      FrequencyType;
  typedef typename ImageType::SpacingValueType FrequencyValueType;
  /** Default constructor. Needed since we provide a cast constructor. */
  FrequencyImageRegionIteratorWithIndex() :
    FrequencyImageRegionConstIteratorWithIndex< TImage >()
  {
  }

  /** Constructor establishes an iterator to walk a particular image and a
   * particular region of that image. */
  FrequencyImageRegionIteratorWithIndex(TImage *ptr, const RegionType & region) :
    FrequencyImageRegionConstIteratorWithIndex< TImage >(ptr, region)
  {
  }

  /** Constructor that can be used to cast from an ImageIterator to an
   * ImageRegionIteratorWithIndex. Many routines return an ImageIterator, but for a
   * particular task, you may want an ImageRegionIteratorWithIndex.  Rather than
   * provide overloaded APIs that return different types of Iterators, itk
   * returns ImageIterators and uses constructors to cast from an
   * ImageIterator to a ImageRegionIteratorWithIndex. */
  FrequencyImageRegionIteratorWithIndex(const ImageIteratorWithIndex< TImage > & it) :
    FrequencyImageRegionConstIteratorWithIndex< TImage >(it)
  {
  }

  /** Set the pixel value */
  void Set(const PixelType & value) const
  { this->m_PixelAccessorFunctor.Set(*( const_cast< InternalPixelType * >( this->m_Position ) ), value); }

  /** Return a reference to the pixel.
   * This method will provide the fastest access to pixel
   * data, but it will NOT support ImageAdaptors. */
  PixelType & Value(void)
  { return *( const_cast< InternalPixelType * >( this->m_Position ) ); }

protected:
  /** The construction from a const iterator is declared protected
      in order to enforce const correctness. */
  FrequencyImageRegionIteratorWithIndex(const FrequencyImageRegionConstIteratorWithIndex< TImage > & it) :
    FrequencyImageRegionConstIteratorWithIndex< TImage >(it)
  {
  }

  Self & operator=(const FrequencyImageRegionConstIteratorWithIndex< TImage > & it)
  {
    this->FrequencyImageRegionConstIteratorWithIndex< TImage >::operator=(it);
    return *this;
  }
};
} // end namespace itk
#endif
