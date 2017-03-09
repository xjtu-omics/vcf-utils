#ifndef EVENT_H
#define EVENT_H

#include <string>

struct Event {
  Event(const std::string& line);

  int chromosomeIndex_;
  int position_;
  std::string refAllele_;
  std::string altAllele_;
};

bool operator<(const Event& leftEvent, const Event& rightEvent);

#endif // EVENT_H
