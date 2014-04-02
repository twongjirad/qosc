#include "EventReweightParameter.hh"

using namespace qosc;

EventReweightParameter::EventReweightParameter( std::string systermname, double mean, double sig ) 
  : BasicParameter( systermname, mean, sig, kMinuitTerm ) {

}

EventReweightParameter::~EventReweightParameter() {
  // who owns the event reweight class?
}

void EventReweightParameter::SetPullValue( double pull ) {
  SetValue( pull );
  std::cout << "EventReweightParameter::SetPullValue(" << pull << ")" << std::endl;
  for ( SysTermIDictIter it=SysTermIDictBegin(); it!=SysTermIDictEnd(); it++ ) {
    // need to pass pull value to eventweight generators
    (*it).second->SetPullValue( pull );
    std::cout << "EventReweightParameter(" << GetName() << "): " << (*it).first << " " << pull << std::endl;
  }
}
