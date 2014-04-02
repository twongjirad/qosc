#include "SysTermI.hh"

using namespace qosc;

SysTermI::SysTermI( std::string name, SysTermI::SysTermType systype ) {
  fType = systype;
  m_name = name;
}

SysTermI::~SysTermI() {
}
