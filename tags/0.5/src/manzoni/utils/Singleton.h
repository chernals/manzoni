#ifndef SINGLETON_H
#define SINGLETON_H

template <class T>
class Singleton
{
	private:
		Singleton();
		~Singleton();
		Singleton(Singleton const&);
		Singleton& operator=(Singleton const&);

	public:
		static T& getInstance()
		{
			static T _instance;
			return _instance;
		}
		static T& get()
		{
			return getInstance();
		}

};

#endif    
