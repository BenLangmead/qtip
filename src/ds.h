#ifndef __qsim__ds__
#define __qsim__ds__

#include <cassert>
#include <string.h>
#include "rnglib.hpp"

template <typename T, int S = 128>
class EList {

public:

	/**
	 * Allocate initial default of S elements.
	 */
	explicit EList() :
		list_(NULL), sz_(S), cur_(0) { }

	/**
	 * Destructor.
	 */
	~EList() { free(); }

	/**
	 * Return number of elements.
	 */
	inline size_t size() const { return cur_; }

	/**
	 * Return number of elements allocated.
	 */
	inline size_t capacity() const { return sz_; }

	/**
	 * Ensure that there is sufficient capacity to expand to include
	 * 'thresh' more elements without having to expand.
	 */
	inline void ensure(size_t thresh) {
		if(list_ == NULL) lazyInit();
		expandCopy(cur_ + thresh);
	}

	/**
	 * Ensure that there is sufficient capacity to include 'newsz' elements.
	 * If there isn't enough capacity right now, expand capacity to exactly
	 * equal 'newsz'.
	 */
	inline void reserveExact(size_t newsz) {
		if(list_ == NULL) lazyInitExact(newsz);
		expandCopyExact(newsz);
	}

	/**
	 * Return true iff there are no elements.
	 */
	inline bool empty() const { return cur_ == 0; }

	/**
	 * Add an element to the back and immediately initialize it via
	 * operator=.
	 */
	void push_back(const T& el) {
		if(list_ == NULL) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		list_[cur_++] = el;
	}

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void expand() {
		if(list_ == NULL) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
	}

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void fill(size_t begin, size_t end, const T& v) {
		assert(begin <= end);
		assert(end <= cur_);
		for(size_t i = begin; i < end; i++) {
			list_[i] = v;
		}
	}

	/**
	 * Add an element to the back.  No intialization is done.
	 */
	void fill(const T& v) {
		for(size_t i = 0; i < cur_; i++) {
			list_[i] = v;
		}
	}

	/**
	 * Set all bits in specified range of elements in list array to 0.
	 */
	void fillZero(size_t begin, size_t end) {
		assert(begin <= end);
		memset(&list_[begin], 0, sizeof(T) * (end-begin));
	}

	/**
	 * Set all bits in the list array to 0.
	 */
	void fillZero() {
		memset(list_, 0, sizeof(T) * cur_);
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resizeNoCopy(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandNoCopy(sz);
		cur_ = sz;
	}

	/**
	 * If size is less than requested size, resize up to at least sz
	 * and set cur_ to requested sz.
	 */
	void resize(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) {
			expandCopy(sz);
		}
		cur_ = sz;
	}

	/**
	 * If size is less than requested size, resize up to exactly sz and set
	 * cur_ to requested sz.
	 */
	void resizeExact(size_t sz) {
		if(sz > 0 && list_ == NULL) lazyInitExact(sz);
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandCopyExact(sz);
		cur_ = sz;
	}

	/**
	 * Erase element at offset idx.
	 */
	void erase(size_t idx) {
		assert(idx < cur_);
		for(size_t i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
	}

	/**
	 * Erase range of elements starting at offset idx and going for len.
	 */
	void erase(size_t idx, size_t len) {
		assert(len >= 0);
		if(len == 0) {
			return;
		}
		assert(idx < cur_);
		for(size_t i = idx; i < cur_-len; i++) {
			list_[i] = list_[i+len];
		}
		cur_ -= len;
	}

	/**
	 * Insert value 'el' at offset 'idx'
	 */
	void insert(const T& el, size_t idx) {
		if(list_ == NULL) lazyInit();
		assert(idx <= cur_);
		if(cur_ == sz_) expandCopy(sz_+1);
		for(size_t i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
	}

	/**
	 * Insert contents of list 'l' at offset 'idx'
	 */
	void insert(const EList<T>& l, size_t idx) {
		if(list_ == NULL) lazyInit();
		assert(idx < cur_);
		if(l.cur_ == 0) return;
		if(cur_ + l.cur_ > sz_) expandCopy(cur_ + l.cur_);
		for(size_t i = cur_ + l.cur_ - 1; i > idx + (l.cur_ - 1); i--) {
			list_[i] = list_[i - l.cur_];
		}
		for(size_t i = 0; i < l.cur_; i++) {
			list_[i+idx] = l.list_[i];
		}
		cur_ += l.cur_;
	}

	/**
	 * Remove an element from the top of the stack.
	 */
	void pop_back() {
		assert(cur_ > 0);
		cur_--;
	}

	/**
	 * Make the stack empty.
	 */
	void clear() {
		cur_ = 0; // re-use stack memory
		// Don't clear heap; re-use it
	}

	/**
	 * Get the element on the top of the stack.
	 */
	inline T& back() {
		assert(cur_ > 0);
		return list_[cur_-1];
	}

	/**
	 * Get the element on the top of the stack, const version.
	 */
	inline const T& back() const {
		assert(cur_ > 0);
		return list_[cur_-1];
	}

	/**
	 * Get the frontmost element (bottom of stack).
	 */
	inline T& front() {
		assert(cur_ > 0);
		return list_[0];
	}

	/**
	 * Get the element on the bottom of the stack, const version.
	 */
	inline const T& front() const { return front(); }

	/**
	 * Return a reference to the ith element.
	 */
	inline T& operator[](size_t i) {
		assert(i < cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline const T& operator[](size_t i) const {
		assert(i < cur_);
		return list_[i];
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline T& get(size_t i) {
		return operator[](i);
	}

	/**
	 * Return a reference to the ith element.
	 */
	inline const T& get(size_t i) const {
		return operator[](i);
	}

	/**
	 * Return a pointer to the beginning of the buffer.
	 */
	T *ptr() { return list_; }

	/**
	 * Return a const pointer to the beginning of the buffer.
	 */
	const T *ptr() const { return list_; }

private:

	/**
	 * Initialize memory for EList.
	 */
	void lazyInit() {
		assert(list_ == NULL);
		list_ = alloc(sz_);
	}

	/**
	 * Initialize exactly the prescribed number of elements for EList.
	 */
	void lazyInitExact(size_t sz) {
		assert(sz > 0);
		assert(list_ == NULL);
		sz_ = sz;
		list_ = alloc(sz);
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	T *alloc(size_t sz) {
		T* tmp = new T[sz];
		assert(tmp != NULL);
		return tmp;
	}

	/**
	 * Allocate a T array of length sz_ and store in list_.  Also,
	 * tally into the global memory tally.
	 */
	void free() {
		if(list_ != NULL) {
			delete[] list_;
			list_ = NULL;
			sz_ = cur_ = 0;
		}
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.  Size
	 * increases quadratically with number of expansions.  Copy old contents
	 * into new buffer using operator=.
	 */
	void expandCopy(size_t thresh) {
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandCopyExact(newsz);
	}

	/**
	 * Expand the list_ buffer until it has exactly 'newsz' elements.  Copy
	 * old contents into new buffer using operator=.
	 */
	void expandCopyExact(size_t newsz) {
		if(newsz <= sz_) return;
		T* tmp = alloc(newsz);
		assert(tmp != NULL);
		size_t cur = cur_;
		if(list_ != NULL) {
			for(size_t i = 0; i < cur_; i++) {
				// Note: operator= is used
				tmp[i] = list_[i];
			}
			memset(list_, 0, sizeof(T) * cur_); // make all ptrs NULL
			free();
		}
		list_ = tmp;
		sz_ = newsz;
		cur_ = cur;
	}

	/**
	 * Expand the list_ buffer until it has at least 'thresh' elements.
	 * Size increases quadratically with number of expansions.  Don't copy old
	 * contents into the new buffer.
	 */
	void expandNoCopy(size_t thresh) {
		assert(list_ != NULL);
		if(thresh <= sz_) return;
		size_t newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandNoCopyExact(newsz);
	}

	/**
	 * Expand the list_ buffer until it has exactly 'newsz' elements.  Don't
	 * copy old contents into the new buffer.
	 */
	void expandNoCopyExact(size_t newsz) {
		assert(list_ != NULL);
		assert(newsz > 0);
		free();
		T* tmp = alloc(newsz);
		assert(tmp != NULL);
		list_ = tmp;
		sz_ = newsz;
		assert(sz_ > 0);
	}

	T *list_;      // list pointer, returned from new[]
	size_t sz_;    // capacity
	size_t cur_;   // occupancy (AKA size)
};

/**
 * Reservoir sampled version of an EList<T>
 */
template<typename T>
class ReservoirSampledEList {
	
public:

	ReservoirSampledEList(size_t k) : k_(k), n_(0), list_() { }

	/**
	 * Possibly add item, using reservoir sampling.
	 */
	void add(const T& t) {
		n_++;
		if(list_.size() < k_) {
			list_.push_back(t);
		} else {
			size_t j = (size_t)(r4_uni_01() * n_);
			assert(j < n_);
			if(j < list_.size()) {
				list_[j] = t;
			}
		}
	}

	/**
	 * Possibly add item, using reservoir sampling.
	 */
	size_t add_part1() {
		n_++;
		if(list_.size() < k_) {
			assert(list_.size() == n_-1);
			list_.expand();
			return list_.size()-1;
		} else {
			return (size_t)(r4_uni_01() * n_);
		}
	}

	/**
	 * Return the number of items added (not all of which were retained by the
	 * sampler).
	 */
	size_t size() const {
		return n_;
	}
	
	/**
	 * Return true iff no items have yet been added.
	 */
	bool empty() const {
		return size() == 0;
	}
	
	/**
	 * Return reservoir size.
	 */
	size_t k() const {
		return k_;
	}
	
	/**
	 * Return const version of list.
	 */
	const EList<T> list() const {
		return list_;
	}

	/**
	 * Return list.
	 */
	EList<T> list() {
		return list_;
	}

protected:
	
	size_t k_;
	size_t n_;
	EList<T> list_;
};


#endif

