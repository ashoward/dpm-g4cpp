#include "TrackStack.hh"

TrackStack& TrackStack::Instance() {
  static TrackStack stack;
  return stack;
}
