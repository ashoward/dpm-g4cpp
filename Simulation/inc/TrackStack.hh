
#ifndef TrackStack_HH
#define TrackStack_HH

//
// A simple (singletone) track-stack to handle both primary and secondary
// particles. Each primary is inserted first, then each history goes till the
// stack becomes empty (i.e. PopIntoThisTrack() returns with -1). During the
// history, seconday tracks can be insterted by using the Insert() method (that
// makes sure that stack is deep enough and calles the Track::Reset() method
// before giving back a reference to the Track.
//

#include "Track.hh"

#include <vector>
#include <iostream>
class TrackStack {
public:
   static TrackStack& Instance();

  // pops a secondary track from the secondary strack and writes to aTrack;
  // and returns with its original index (empty stack if returns -1)
  int PopIntoThisTrack(Track& track) {
    // return -1 if the secondary stack is empty
    if (fCurIndx<0) {
      return -1;
    }
    // compy the next avaiable seconday track to the primary
    fSecondaries[fCurIndx].Copy(track);
    // return with the currently used secondary index and decrease
    return fCurIndx--;
  }

  // returns a reference to a secondary track (that is reset before)
  // NOTE: the stack is resized when needed
  Track&  Insert() {
    // make sure that the size if fine
    ++fCurIndx;
    if (fCurIndx==fSize) {
      fSize *= 2;
      fSecondaries.resize(fSize);
    }
    // retrun a eference to the next avaiable secondary track
    fSecondaries[fCurIndx].Reset();
    return fSecondaries[fCurIndx];
  }

private:
  TrackStack() {
    fSize    = 16;
    fCurIndx = -1;
    fSecondaries.resize(fSize);
  }

private:
  int fSize;                       // #tracks in the stack
  int fCurIndx;                    // the (max) used seconday index on the stack
  std::vector<Track> fSecondaries; // secondary tracks
};

#endif // TrackStack_HH
