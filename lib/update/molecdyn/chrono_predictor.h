#ifndef chrono_predictor_h
#define chrono_predictor_h

namespace Chroma {

  // Abstract interface for a Chronological Solution predictor
  template<typename T>
  class AbsChronologicalPredictor {
  public:

    // Set psi to be the next initial guess
    virtual void operator()(T& psi) const = 0;

    // Present new vector for use in future chronological
    // Predictors
    virtual void newVector(const T& psi) = 0;
  };


  // A Special Chronological inverter that is not chronological
  // at all but always provides a zero guess
  template<typename T> 
  class ZeroGuess : public AbsChronologicalPredictor<T> {
  public:
    // Set T to the next initial guess
    void operator()(T& psi) const {
      // This works for multi1d-s too right ?
      psi = zero;
    };

    // This is just for form... Doesnt do anything
    void newVector(const T& psi) {}
  }
}; // End namespace
#endif
