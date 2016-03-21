#ifndef CACHE_H
#define CACHE_H

#include <iostream>
#include <map>
#include <stdlib.h>
using namespace std;

class CacheNode{
	public:
	CacheNode(int ind, double* val){
		index = ind;
		value = val;
		prev = NULL;
		next = NULL;
	}
	void dumpInfo(ostream& out){
		out << "prev=" << prev << ", next=" << next << ", index=" << index ;
	}
	CacheNode* prev;
	CacheNode* next;
	int index;
	double* value;
};

class ArrayCache{
	
	public:
	ArrayCache(int _size_limit){
		size_limit = _size_limit;
		head = NULL;
		tail = NULL;
	}
	
	//return false if not inserted due to duplicate
	bool insert(int ind, double* arr){
		//check if exist, do nothing
		map<int,CacheNode*>::iterator it;
	       	if( (it=cache_map.find(ind)) != cache_map.end() )
			return false;
		
		//create CacheNode for the element, put at the head
		CacheNode* node = new CacheNode(ind, arr);
		addToHead(node);
		
		//put the node into map
		cache_map.insert( make_pair(ind, node) );
		
		//check if reach limit. If yes, remove tail node
		if( cache_map.size() > size_limit ){
			removeTail();
		}
		
		return true;
	}
		
	//return NULL if not in the cache
	double* get(int ind){
		map<int,CacheNode*>::iterator it;
		if( (it=cache_map.find(ind)) != cache_map.end() ){
			CacheNode* node = it->second;
			moveToHead(node);
			return node->value;
		}else{
			return NULL;
		}
	}

	void dumpInfo(ostream& out){
		out << "size_limit=" << size_limit << endl;
		//out << "head addr=" << head << endl;
		//out << "tail addr=" << tail << endl;
		out << "map.size()=" << cache_map.size() << endl;
		/*out << "cache_map:" << endl;
		for(map<int,CacheNode*>::iterator it=cache_map.begin(); it!=cache_map.end(); it++){
			out << it->first << ", address=" << it->second << ", ";
			it->second->dumpInfo(out);
			out << endl;
		}*/
	}
	
	int size(){
		return cache_map.size();
	}

	private:
	
	void addToHead(CacheNode* node){
		
		if( head!=NULL ){
			CacheNode* tmp = head;
			node->next = tmp;
			tmp->prev = node;
			head = node;
		}else{
			head = node;
			tail = node;
		}
	}

	void moveToHead(CacheNode* node){
		
		//detach node from the current position
		//dumpInfo(cerr);
		detachFromList(node);
		//attach this node to head
		//dumpInfo(cerr);
		addToHead(node);
	}

	void removeTail(){
		if( tail==NULL )
			return ;
		
		CacheNode* node = tail;
		//detach this node from the linkedlist
		detachFromList(node);
		//remove this node from map
		int ind = node->index;
		cache_map.erase(ind);
		//remove this node (and its array) from the memory
		delete[] node->value;
		delete node;
	}
	
	void detachFromList(CacheNode* node){
		
		//cerr << "detach node " << node << endl;

		if( node==tail )
			tail = node->prev;
		if( node==head )
			head = node->next;
		
		if( node->prev != NULL ){ // if the node has prev
			node->prev->next = node->next;
		}
		if(node->next != NULL){ //if the node has next
			node->next->prev = node->prev;
		}
			
		node->prev = NULL;
		node->next = NULL;
	}
	
	int size_limit;
	map<int,CacheNode*> cache_map;
	CacheNode* head;
	CacheNode* tail;
};

#endif
