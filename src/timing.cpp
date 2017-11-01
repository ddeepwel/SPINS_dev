#include "timing.hpp"
#include <string>
#include <vector>
#include <mpi.h>
#include <stdio.h>

#ifdef TIMING_ENABLE
// Structure for information related to a single timing invocation
struct timing_element; 

struct timing_element {
   std::string name; //name of the current element
   int count; // number of times this element has been pushed
   double tic_time; // MPI_Wtime when this was pushed to the stack
   double total_time; // Total accumulated runtime so far
   std::vector<timing_element *> children; // List of sub-elements pushed thus far, at any time
   timing_element *parent; // Pointer to the parent element
};

std::vector<timing_element *> root(0); // Elements created at root level

timing_element *current = 0; // Pointer to current element

/* Add a new element to the timing stack.  Re-use an existing element
   of this name if one already exists, otherwise create a new structure
   and append it to the current-children list */
void timing_push(const char * name){
   //fprintf(stderr,"Pushing new element...\n");
   std::vector<timing_element *> *child_list; // list to check for duplicates
   if (current == 0) {
      child_list = &root;
   } else {
      child_list = &current->children;
   }

   for (int i = 0; i < child_list->size(); i++) {
     // fprintf(stderr,"Considering element %d, name %s\n",i,(*child_list)[i]->name.c_str());

      if (!strcmp((*child_list)[i]->name.c_str(),name)) {
         // We have a match, so move into this element
         timing_element *child = (*child_list)[i];
         child->tic_time = MPI_Wtime();
         child->count += 1;
         current = child;
         return; // done
      }
   }

   // no match found, so construct a new element
   /*fprintf(stderr,"Constructing new element %s ",name);
   if (current == 0) {
      fprintf(stderr,"as a child of root\n");
   } else {
      fprintf(stderr,"as a child of element %s\n",current->name.c_str());
   }*/

   timing_element *new_el = new timing_element;
   new_el->name = name;
   new_el->tic_time = MPI_Wtime();
   new_el->total_time = 0;
   new_el->count = 1;
   new_el->parent = current;

   child_list->push_back(new_el);
   current = new_el;
   return;
}

void timing_pop(){
   if (current == 0) {
      fprintf(stderr,"ERROR: Popping at the top of the timing stack!\n");
      return;
   } 
   double now = MPI_Wtime();
   //fprintf(stderr,"Popping timing element %s\n",current->name.c_str());
   current->total_time += now - current->tic_time;
   current = current->parent;
}

void print_element(timing_element * el, int indent) {
   if (el->total_time < 1) { // output in ms
      fprintf(stderr,"%*s%s -- %d -- %.2fms\n",indent,"",el->name.c_str(),el->count,el->total_time*1000);
   } else { // output in seconds
      fprintf(stderr,"%*s%s -- %d -- %.3fs\n",indent,"",el->name.c_str(),el->count,el->total_time);
   }
}

void timing_stack_recurse(timing_element * el, int indent) {
   print_element(el, indent);
   for (int i = 0; i < el->children.size(); i++) {
      timing_stack_recurse(el->children[i], indent + 1);
   }
}

void timing_stack_report(){
   for (int i = 0; i < root.size(); i++) {
      timing_stack_recurse(root[i],0);
   }
}

#else

void timing_stack_report() {
   fprintf(stderr,"Timing was not enabled for this build of SPINS.\n");
}

#endif
