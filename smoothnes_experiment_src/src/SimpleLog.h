#include "Tools.h"
#include <fstream> // std::ofstream

#ifndef SIMPLELOG_H_
#define SIMPLELOG_H_

class SimpleLog
{
  public:
	string filename;

	string getCurrentDateTime(string s)
	{
		time_t now = time(0);
		struct tm tstruct;
		char buf[80];
		tstruct = *localtime(&now);
		if (s == "now")
			strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
		else if (s == "date")
			strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);
		return string(buf);
	};

	template <class T>
	void write(T message, bool newline = true)
	{

		ofstream myfile;
		myfile.open(filename, ios::out | ios::app);
		myfile << toString(message);
		if (newline)
		{
			myfile << '\n';
		}
		myfile.close();
	}

	template <class T>
	void write_variable(string variable_name, T variable, bool newline = true)
	{

		ofstream myfile;
		myfile.open(filename, ios::out | ios::app);
		myfile << variable_name << ": ";
		myfile << toString(variable);
		if (newline)
		{
			myfile << '\n';
		}
		myfile.close();
	}

	SimpleLog(string filename, bool write_start = true)
	{

		this->filename = filename;
		if (write_start)
		{
			this->write("-----------------------------------------");
			this->write(this->getCurrentDateTime("now"));
		}
	}
	~SimpleLog()
	{
	}

	template <class T>
	void write_vector(T *vector, int len, string vector_name = "0")
	{
		if (vector_name == "0")
		{
			for (int i = 0; i < len; i++)
			{
				if (i < len - 1)
				{
					this->write(toString(vector[i]) + " ", false);
				}
				else
				{
					this->write(toString(vector[i]), true);
				}
			}
		}
		else
		{
			this->write(vector_name + ": ", false);
			this->write_vector(vector, len, "0");
		}
	}
};

#endif
