/**
* @file   blist.h
* @author Jed Brown <jed@59A2.org>
* @date   Sun Aug 24 20:36:47 2008
* 
* @brief A simple way of managing memory in chunks.
*
* It is just a simple linked list of memory chunks.  Ask for more memory when you need it.  Free the whole thing at
* once.  It's better to use a few big chunks than lots of small ones.
* 
*/
#ifndef _BLIST_H
#define _BLIST_H

#include "dohptype.h"

typedef struct m_dBufferList *dBufferList;
struct m_dBufferList {
  void *data;
  dBufferList next;
};

static dErr dBufferListFree(dBufferList *head)
{
  dBufferList cur,next;
  dErr err;

  dFunctionBegin;
  dValidPointer(head,1);        /* confirm that head is a valid pointer, not that there is actually a link */
  cur = *head;
  while (cur) {
    next = cur->next;
    err = dFree(cur->data);dCHK(err);
    err = dFree(cur);dCHK(err);
    cur = next;
  }
  dFunctionReturn(0);
}

static dErr dBufferListMalloc(dBufferList *head,size_t bytes,void **data)
{
  dBufferList cur,*addr;
  dErr err;

  dFunctionBegin;
  dValidPointer(head,1);
  dValidPointer(data,3);
  if (!bytes) dFunctionReturn(0);
  addr = head;
  cur = *head;
  while (cur && cur->data) { /* Walk the list until we find the end or a node with no data. */
    addr = &cur->next;
    cur = cur->next;
  }
  if (!cur) {             /* We found the end of the list, allocate a new node in \c addr */
    err = dNew(struct m_dBufferList,addr);dCHK(err);
    cur = *addr;
  }
  /* Either way, we are now at a node with empty data. */
  err = dMalloc(bytes,&cur->data);dCHK(err);
  *data = cur->data;
  dFunctionReturn(0);
}

#endif
