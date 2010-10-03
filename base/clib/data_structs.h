/********************************************************************
 *
 *  data_struct.h
 *
 *  Template classes of some basic data structure.
 *
 *  Created by Dahua Lin, on Oct 03, 2010
 *
 ********************************************************************/


namespace smi
{
    
template<typename T>
struct RefMemory
{
    int n;
    T *base;
    
    RefMemory(int n_, T *p) : n(n_), base(p) { }
};
    

template<typename T>
RefMemory<T> refmem(int n, T *p)
{
    return RefMemory<T>(n, p);
}

    
template<typename T>
class SeqList
{
public:
    explicit SeqList(int n) 
    : m_capa(n), m_n(0), m_data(new T[n]), m_own(true)
    {
    }
    
    explicit SeqList(RefMemory<T> r)
    : m_capa(r.n), m_n(0), m_data(r.base), m_own(false)
    {
    }
    
    ~SeqList()
    {
        if (m_own)
            delete[] m_data;
    }    
    
    int capacity() const
    {
        return m_capa;
    }
    
    int size() const
    {
        return m_n;
    }
    
    bool empty() const
    {
        return m_n == 0;
    }
    
    
    const T *data() const
    {
        return m_data;
    }
        
    T *data()
    {
        return m_data;
    }
        
    const T& operator[] (int i) const
    {
        return m_data[i];
    }
    
    T& operator[] (int i) 
    {
        return m_data[i];
    }
    
    void add(const T& e) 
    {
        m_data[m_n++] = e;
    }    
    
private:
    int m_capa;
    int m_n;
    T *m_data;
    bool m_own;
    
}; // end class SeqList
    


template<typename T>
class Stack
{
public:
    explicit Stack(int n) 
    : m_capa(n), m_n(0), m_data(new T[n]), m_own(true)
    {
    }
    
    explicit Stack(RefMemory<T> r)
    : m_capa(r.n), m_n(0), m_data(r.base), m_own(false)
    {
    }
    
    ~Stack()
    {
        if (m_own)
            delete[] m_data;
    }    
    
    int capacity() const
    {
        return m_capa;
    }
    
    int size() const
    {
        return m_n;
    }
    
    bool empty() const
    {
        return m_n == 0;
    }
    
    void push(const T& e) 
    {
        m_data[m_n++] = e;
    }
    
    void pop()
    {
        -- m_n;
    }
    
    const T& top() const
    {
        return m_data[m_n];
    }
    
    T& top() 
    {
        return m_data[m_n];
    }
    
private:
    int m_capa;
    int m_n;
    T *m_data;
    bool m_own;
    
}; // end class Stack



  

template<typename T>
class Queue
{
public:
    explicit Queue(int n)
    : m_capa(n), m_ifront(0), m_iback(0), m_data(new T[n]), m_own(true)
    {
    }
    
    explicit Queue(RefMemory<T> r)
    : m_capa(r.n), m_ifront(0), m_iback(0), m_data(r.base), m_own(false)
    {
    }
    
    ~Queue()
    {
        if (m_own)
            delete[] m_data;
    }    
    
    int capacity() const
    {
        return m_capa;
    }
    
    int size() const
    {
        return m_iback - m_ifront;
    }
    
    bool empty() const
    {
        return m_ifront == m_iback;
    }
    
    void enqueue(const T& e)
    {
        m_data[m_iback++] = e;
    }
    
    void dequeue()
    {
        m_ifront++;
    }
    
    const T& front() const
    {
        return m_data[m_ifront];
    }
    
    T& front()
    {
        return m_data[m_ifront];
    }            
    
private:
    int m_capa;
    int m_ifront;
    int m_iback;
    T *m_data;
    bool m_own;
    
}; // end class Queue

    
}




