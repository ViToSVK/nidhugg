#ifndef NDEBUG

#include <iostream>
#include <sstream>

extern const bool DEBUG;

static int err_level = 0;

static void indent(bool group = DEBUG) {
  if (!group) return;
  for (int i = 0; i < err_level; i++) {
    std::cerr << "  ";
  }
}

static void err_msg(const std::string& str, bool group = DEBUG) {
  if (!group) return;
  indent();
  std::cerr << str << "\n";
}

static void start_err(const std::string& str, bool group = DEBUG) {
  if (!group) return;
  err_msg(str);
  err_level += 1;
}

static void end_err(const std::string& str = "", bool group = DEBUG) {
  if (!group) return;
  assert(err_level > 0);
  err_level -= 1;
  std::stringstream ss;
  ss << "Done. (" << str << ")";
  err_msg(ss.str());
}

#else

/* What a hack! When NDEBUG is defined, all of the debug-functions are
 * interpreted as macros that expand into nothing. */

#define err_msg(...);
#define start_err(...);
#define end_err(...);

#endif
