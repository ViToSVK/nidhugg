/* Copyright (C) 2016-2017 Marek Chalupa
 * Copyright (C) 2017-2019 Viktor Toman
 *
 * This file is part of Nidhugg.
 *
 * Nidhugg is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Nidhugg is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

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
  std::stringstream ss;
  ss << "Done. (" << str << ")";
  err_msg(ss.str());
  err_level -= 1;
}
