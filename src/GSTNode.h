/*
 * GSTNode.h
 *
 *  Created on: Apr 19, 2011
 *      Author: mgymrek
 */

#ifndef GSTNODE_H_
#define GSTNODE_H_
#include <list>

using namespace std;

/*
 * (self.ids[ident],index,len(seq)
 */
struct Alignment{
	int ident;
	int index;
	int length;

};

/*
 * {"seq":node.labelstr(),"alignments":[]}
 */
struct GSTNode {
	list<Alignment> alignments;
};

#endif /* GSTNODE_H_ */
