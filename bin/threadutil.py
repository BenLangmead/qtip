import threading
try:
    from Queue import Queue, Empty, Full
except ImportError:
    from queue import Queue, Empty, Full  # python 3.x

class AsyncWriter(threading.Thread):
    
    ''' Thread that takes items off a queue and writes them to a named file,
        which could be a FIFO '''
    
    def __init__(self, fn, q):
        threading.Thread.__init__(self)
        self.fn = fn
        self.q = q
    
    def run(self):
        ''' Write reads to the FIFO.  When None appears on the queue, we're
            done receiving reads and we close the FIFO filehandle. '''
        i = 0
        fh = open(self.fn, 'w')
        while True:
            item = self.q.get()
            if item is None:
                # None signals the final write
                self.q.task_done()
                fh.close()
                break
            fh.write(item)
            self.q.task_done()
            i += 1
        fh.close() # close FIFO filehandle

def q2fh(q, fh, timeout=0.2):
    ''' Get strings from 'queue' and write them to filehandle 'inp'.
        When we dequeue None, we're done. '''
    while True:
        try:
            i = q.get(block=True, timeout=timeout)
        except Empty:
            continue
        if i is None:
            break
        fh.write(i)
    fh.close()

def fh2q(fh, q, timeout=0.2):
    ''' Get strings from 'out' filehandle and add them to 'queue'.  We
        enqueue None last, to indicate no more inputs. '''
    for ln in fh:
        while True:
            try:
                q.put(ln, block=True, timeout=timeout)
                break
            except Full:
                continue
    q.put(None)
    fh.close()
