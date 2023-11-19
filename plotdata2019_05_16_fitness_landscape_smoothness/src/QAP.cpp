/*
 * hamming_functions.cpp
 *
 *  Created on: Jun 6, 2018
 */
#include "QAP.h"
#include <assert.h>
#include "Tools.h"
#include <assert.h>
#include <iostream>
QAP::QAP()
{
}

QAP::~QAP()
{
	for (int i = 0; i < n; i++)
	{
		delete[] m_distance_matrix[i];
		delete[] m_flow_matrix[i];
	}
	delete[] m_flow_matrix;
	delete[] m_distance_matrix;
	delete[] m_aux;
}

void QAP::read(string filename)
{

	char line[5096]; // variable for input value
	ifstream indata;

	indata.open(filename.c_str(), ios::in);

	assert(indata.is_open());

	int num = 0;
	//int row=0;
	//int col=0;
	while (!indata.eof())
	{
		//LEER LA LINEA DEL FICHERO
		indata.getline(line, 5096);
		stringstream ss;
		string sline;
		ss << line;
		ss >> sline;
		if (num == 0)
		{
			//OBTENER EL TAMAÑO DEL PROBLEMA
			n = atoi(sline.c_str());

			assert(n > 1);

			m_distance_matrix = new int *[n];
			m_flow_matrix = new int *[n];
			for (int i = 0; i < n; i++)
			{
				m_distance_matrix[i] = new int[n];
				m_flow_matrix[i] = new int[n];
			}
			//row=0;
			//col=0;
		}
		else if (1 <= num && num <= n)
		{
			//LOAD DISTANCE MATRIX
			char *pch;
			pch = strtok(line, " ");
			int distance = atoi(pch);
			m_distance_matrix[num - 1][0] = distance;
			for (int i = 1; i < n; i++)
			{
				pch = strtok(NULL, " ,.");
				distance = atoi(pch);
				m_distance_matrix[num - 1][i] = distance;
			}
		}
		else if (num > n && num <= (2 * n))
		{
			//LOAD FLOW MATRIX
			char *pch;
			pch = strtok(line, " ");
			int weight = atoi(pch);
			m_flow_matrix[num - n - 1][0] = weight;
			for (int i = 1; i < n; i++)
			{
				pch = strtok(NULL, " ,.");
				weight = atoi(pch);
				m_flow_matrix[num - n - 1][i] = weight;
			}
		}

		else
		{
			break;
		}

		/*
        //LOAD DISTANCE MATRIX
        else{
            int distance;
            if (row>=m_size){
                //flow_matrix
                char * pch;
                pch = strtok (line," ");
                while (pch != NULL)
                {
                    distance=atoi(pch);
                    m_flow_matrix[row-m_size][col]=distance;
                    //cout<<"flow: "<<distance<<" row: "<<row<<" col: "<<col<<endl;
                    col++;
                    if (col==m_size){
                        col=0;
                        row++;
                    }
                    pch = strtok (NULL, " ,.");
                }
            }
            else{
                //distance_matrix
                char * pch;
                pch = strtok (line," ");
                while (pch != NULL)
                {
                    distance=atoi(pch);
                 //   cout<<"dist: "<<distance<<" row: "<<row<<" col: "<<col<<endl;
                    m_distance_matrix[row][col]=distance;
                    col++;
                    if (col==m_size){
                        col=0;
                        row++;
                    }
                    pch = strtok (NULL, " ,.");
                }
            }
        }*/

		num++;
	}
	//PrintMatrix(m_distance_matrix, m_size, m_size, "");
	//PrintMatrix(m_flow_matrix, m_size, m_size, "");
	//exit(1);
	indata.close();
	m_aux = new int[n];
}

void QAP::read(void)
{
	string s_line;
	char line[5096];
	int read_line_index = 0;

	while (getline(cin, s_line))
	{
		strcpy(line, s_line.c_str());

		if (read_line_index == 0)
		{
			//OBTENER EL TAMAÑO DEL PROBLEMA
			n = atoi(line);
			cout << "Problem size: " << n << endl;

			assert(n > 1);

			m_distance_matrix = new int *[n];
			m_flow_matrix = new int *[n];
			for (int i = 0; i < n; i++)
			{
				m_distance_matrix[i] = new int[n];
				m_flow_matrix[i] = new int[n];
			}
			//row=0;
			//col=0;
		}
		else if (1 <= read_line_index && read_line_index <= n)
		{
			//LOAD DISTANCE MATRIX
			char *pch;
			pch = strtok(line, " ");
			int distance = atoi(pch);
			m_distance_matrix[read_line_index - 1][0] = distance;
			for (int i = 1; i < n; i++)
			{
				pch = strtok(NULL, " ,.");
				distance = atoi(pch);
				m_distance_matrix[read_line_index - 1][i] = distance;
			}
		}
		else if (read_line_index > n && read_line_index <= (2 * n))
		{
			char *pch;
			pch = strtok(line, " ");
			int weight = atoi(pch);
			m_flow_matrix[read_line_index - n - 1][0] = weight;
			for (int i = 1; i < n; i++)
			{
				pch = strtok(NULL, " ,.");
				weight = atoi(pch);
				m_flow_matrix[read_line_index - n - 1][i] = weight;
			}
		}

		else
		{
			break;
		}

		read_line_index++;
	}
	m_aux = new int[n];
}

int QAP::evaluate(int *sigma)
{
	int fitness = 0;
	int distAB, flowAB, i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			//Hau ondo dago
			distAB = m_distance_matrix[i][j];
			flowAB = m_flow_matrix[sigma[i] - 1][sigma[j] - 1];
			fitness = fitness + (distAB * flowAB);
			//cout << fitness << endl;
			//assert(fitness >= 0);
		}
	}

	return -fitness;
}

int QAP::update_fitness_on_swap(int *sigma_1, int *sigma_2)
{
	int new_fitness_delta = 0;
	int pos1, pos2;

	for (int i = 0; i < n; i++)
	{
		if (sigma_1[i] != sigma_2[i])
		{
			pos1 = i;
			for (int j = i + 1; j < n; j++)
			{
				if (sigma_1[j] != sigma_2[j])
				{
					pos2 = j;
					break;
				}
			}
			break;
		}

		if (i == n)
		{
			//reaching this part of the code means that sigma_1 and sigma_2 are equal, so they should have the same fitness
			return 0;
		}
	}

	// cout << "--QAP->swap_fitness_START--" << endl;
	// cout << "Pos1, pos2:" << pos1 << ", " << pos2 << endl;
	// cout << evaluate(sigma_1) << " - ";
	// PrintPythonArray(sigma_1, n);
	// cout << evaluate(sigma_2) << " - ";
	// PrintPythonArray(sigma_2, n);

	for (int i = 0; i < n; i++)
	{

		new_fitness_delta += m_distance_matrix[i][pos1] * (m_flow_matrix[sigma_2[i] - 1][sigma_2[pos1] - 1] - m_flow_matrix[sigma_1[i] - 1][sigma_1[pos1] - 1]);

		new_fitness_delta += m_distance_matrix[pos1][i] * (m_flow_matrix[sigma_2[pos1] - 1][sigma_2[i] - 1] - m_flow_matrix[sigma_1[pos1] - 1][sigma_1[i] - 1]);

		if (i == pos1)
		{
			continue;
		}

		new_fitness_delta += m_distance_matrix[i][pos2] * (m_flow_matrix[sigma_2[i] - 1][sigma_2[pos2] - 1] - m_flow_matrix[sigma_1[i] - 1][sigma_1[pos2] - 1]);

		new_fitness_delta += m_distance_matrix[pos2][i] * (m_flow_matrix[sigma_2[pos2] - 1][sigma_2[i] - 1] - m_flow_matrix[sigma_1[pos2] - 1][sigma_1[i] - 1]);
	}

	// cout << "new_fitness_delta: " << new_fitness_delta << endl;
	// cout << "fatf_: " << new_fitness_delta + evaluate(sigma_1) << endl;

	// cout << "--QAP->swap_fitness_END--" << endl;

	return -new_fitness_delta;
}
