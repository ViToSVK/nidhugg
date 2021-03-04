/* Copyright (C) 2014-2017 Carl Leonardsson
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

#include "CPid.h"

#include <cassert>
#include <sstream>

CPid::CPid() : aux_idx(-1) {
  _hash = compute_hash();
}

CPid::CPid(const std::vector<int> &pvec){
  assert(pvec.size() > 0);
  assert(pvec[0] == 0);
  auto b = pvec.begin();
  ++b;
  proc_seq = std::vector<int>(b,pvec.end());
  aux_idx = -1;
  _hash = compute_hash();
}

CPid::CPid(const std::initializer_list<int> &il){
  auto b = il.begin();
  assert(b != il.end());
  assert(*b == 0);
  ++b;
  proc_seq = std::vector<int>(b,il.end());
  aux_idx = -1;
  _hash = compute_hash();
}

CPid::CPid(const std::vector<int> &pvec, int i){
  assert(pvec.size() > 0);
  assert(pvec[0] == 0);
  assert(i >= 0);
  auto b = pvec.begin();
  ++b;
  proc_seq = std::vector<int>(b,pvec.end());
  aux_idx = i;
  _hash = compute_hash();
}

CPid CPid::spawn(int pn1) const{
  assert(!is_auxiliary());
  assert(0 <= pn1);
  CPid c = *this;
  c.proc_seq.push_back(pn1);
  c._hash = c.compute_hash();
  return c;
}

CPid CPid::aux(int i) const{
  assert(!is_auxiliary());
  assert(0 <= i);
  CPid c = *this;
  c.aux_idx = i;
  c._hash = c.compute_hash();
  return c;
}

bool CPid::is_auxiliary() const{
  return aux_idx >= 0;
}

std::string CPid::to_string() const{
  std::stringstream ss;
  ss << "<0";
  for(unsigned i = 0; i < proc_seq.size(); ++i){
    ss << "." << proc_seq[i];
  }
  if(is_auxiliary()){
    ss << "/" << aux_idx;
  }
  ss << ">";
  return ss.str();
}

CPid CPid::parent() const{
  assert(has_parent());
  CPid cp(*this);
  if(is_auxiliary()){
    cp.aux_idx = -1;
  }else{
    cp.proc_seq.pop_back();
  }
  cp._hash = cp.compute_hash();
  return cp;
}

bool CPid::has_parent() const{
  return proc_seq.size() || is_auxiliary();
}

std::size_t CPid::get_hash() const{
  assert(compute_hash() == _hash);
  return _hash;
}

std::size_t CPid::compute_hash() const{
  if (proc_seq.empty())
    return (std::size_t) (aux_idx + 1);
  std::size_t res = (proc_seq.size() % 5) * 100;
  if (aux_idx < 0) {
    assert(aux_idx == -1);
    res += 499;
  } else
    res += (aux_idx % 500);
  assert(res <= 999);
  res *= 100;
  std::size_t max = 0;
  for (unsigned i = 1; i < proc_seq.size(); ++i)
    if (proc_seq[i] > 0 && (std::size_t) proc_seq[i] > max)
      max = (std::size_t) proc_seq[i];
  res += ((proc_seq[0] + max) % 100);
  assert(res <= 99999);
  return res;
}

int CPid::compare(const CPid &c) const{
  std::size_t h = get_hash();
  std::size_t ch = c.get_hash();
  if(h < ch) return -1;
  if(h > ch) return 1;

  auto bad_hash = [&]()
  {
    /*
    llvm::errs() << "CPID hash clash: ";
    llvm::errs() << to_string() << " :: " << get_hash() << " --- ";
    llvm::errs() << c.to_string() << " :: " << c.get_hash() << "\n";
    assert(false && "CPID hash clash");
    */
  };

  unsigned i = 0;
  while(i < proc_seq.size() && i < c.proc_seq.size()){
    if(proc_seq[i] < c.proc_seq[i]) { bad_hash(); return -1; }
    if(proc_seq[i] > c.proc_seq[i]) { bad_hash(); return 1; }
    ++i;
  }
  if(i < c.proc_seq.size()) { bad_hash(); return -1; }
  if(i < proc_seq.size()) { bad_hash(); return 1; }
  if(aux_idx == c.aux_idx || (!is_auxiliary() && !c.is_auxiliary())){
    return 0;
  }
  if(aux_idx < c.aux_idx) { bad_hash(); return -1; }
  bad_hash(); return 1;
}

int CPid::get_aux_index() const{
  assert(is_auxiliary());
  return aux_idx;
}

std::vector<int> CPid::get_proc_seq() const{
  return std::vector<int>(proc_seq);
}

CPidSystem::CPidSystem(){
  real_children.push_back({});
  aux_children.push_back({});
  parent.push_back(-1);
  cpids.push_back(CPid());
  identifiers[CPid()] = 0;
}

CPid CPidSystem::dry_spawn(const CPid &c){
  int c_id = identifiers[c];
  return c.spawn(real_children[c_id].size());
}

CPid CPidSystem::spawn(const CPid &c){
  int c_id = identifiers[c];
  CPid c2 = c.spawn(real_children[c_id].size());
  int c2_id = cpids.size();
  real_children[c_id].push_back(c2_id);
  real_children.push_back({});
  aux_children.push_back({});
  parent.push_back(c_id);
  cpids.push_back(c2);
  identifiers[c2] = c2_id;
  return c2;
}

CPid CPidSystem::new_aux(const CPid &c){
  int c_id = identifiers[c];
  CPid c2 = c.aux(aux_children[c_id].size());
  int c2_id = cpids.size();
  aux_children[c_id].push_back(c2_id);
  real_children.push_back({});
  aux_children.push_back({});
  parent.push_back(c_id);
  cpids.push_back(c2);
  identifiers[c2] = c2_id;
  return c2;
}
