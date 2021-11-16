
#include <iostream>
#include <mpi.h>


// Задание 1. Hello world из всех процессов.
void task1() 
{
	int rank, size, recv;
	MPI_Status st;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (!rank) 
	{
		printf("Hello from process %d out of %d!\n", rank, size);
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&recv, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
				MPI_COMM_WORLD, &st);
			printf("Hello from process %d out of %d!\n", recv, size);			
		}
	}
	else
	{
		MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
}

// Задание 2. Максимум массива
void task2()
{
	const int N = 20;
	int a[N], rank, size;
	int max = 0, pmax = 0, maxTrue = 0;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (!rank) 
	{
		for (int i = 0; i < N; i++)
		{
			a[i] = rand() % 100;
			printf("%d ", a[i]);
			if (maxTrue < a[i]) maxTrue = a[i];
		}
		printf("\n");
	}

	int portion = N / size;
	int i1 = portion * rank;
	int i2 = portion * (rank + 1);
	if (rank == size - 1)
		i2 = N;

	MPI_Bcast(a, N, MPI_INT, 0, MPI_COMM_WORLD);

	pmax = 0;
	for (int i = i1; i < i2; i++) {
		if (pmax < a[i]) {
			pmax = a[i];
		}
	};

	MPI_Reduce(&pmax, &max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (!rank) {
		printf("MaxTrue = %d  max = %d \n", maxTrue, max);
	}
	
	MPI_Finalize();

}

// Задание 3. Пи методом монте карло
void task3()
{
	int rank, size;
	double x, y;
	long int localCount, count, i, portion, N = 100000;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	count = 0;
	localCount = 0;
	srand((unsigned)(time(0)));
	portion = N / size;
	for (i = 1; i <= portion; i++) {
		x = ((double)rand()) / ((double)RAND_MAX);
		y = ((double)rand()) / ((double)RAND_MAX);
		if (((x * x) + (y * y)) <= 1) localCount++;
	}
	MPI_Reduce(&localCount, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	double est;
	if (!rank){
		est = ((double)count * 4) / (double)N;
		printf("Pi: %f", est);
	}
	MPI_Finalize();
}

//Задание 4
void task4()
{
	const int N = 20;
	int size, rank, i, a[N];

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (!rank) {
		for (i = 0; i < N; i++) {
			a[i] = rand() % 41 - 20;
			printf("%d ", a[i]);
		}
		printf("\n");
	}

	int* len = new int[size];
	int* ind = new int[size];
	int rest = N;
	int portion = rest / size;

	ind[0] = 0;
	len[0] = portion;
	for (int i = 1; i < size; i++)
	{
		rest -= portion;
		portion = rest / (size - i);
		len[i] = portion;
		ind[i] = ind[i - 1] + len[i - 1];
	}

	int procSize = len[rank];
	int* procA = new int[procSize];

	MPI_Scatterv(a, len, ind, MPI_INT, procA, procSize, MPI_INT, 0, MPI_COMM_WORLD);

	int procSum = 0, sum = 0, procNum = 0, num = 0;
	
	for (int i = 0; i < procSize; i++)
		if (procA[i] > 0)
		{
			procSum += procA[i];
			procNum++;
		}

	MPI_Reduce(&procSum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&procNum, &num, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (!rank) {
		double res = 0.0;
		res = (1.0 * sum) / (1.0 * num);
		printf("Answer = %f \n", res);
	}
	MPI_Finalize();

}


void task5() {
	int size, rank, i;
	const int N = 5;
	int a[N], b[N];
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (!rank) {
		for (i = 0; i < N; i++) {
			a[i] = rand() % 21 - 10;
			b[i] = rand() % 31 - 10;
			printf("%d %d \n", a[i], b[i]);
		}
	}
	int* len = new int[size];
	int* ind = new int[size];
	int rest = N;
	int portion = rest / size;

	ind[0] = 0;
	len[0] = portion;
	for (int i = 1; i < size; i++)
	{
		rest -= portion;
		portion = rest / (size - i);
		len[i] = portion;
		ind[i] = ind[i - 1] + len[i - 1];
	}

	int procSize = len[rank];
	int* procA = new int[procSize];
	int* procB = new int[procSize];


	MPI_Scatterv(a, len, ind, MPI_INT, procA, procSize, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(b, len, ind, MPI_INT, procB, procSize, MPI_INT, 0, MPI_COMM_WORLD);

	int procSum = 0, sum = 0;

	for (int i = 0; i < procSize; i++) {
		procSum += procA[i] * procB[i];
	}

	MPI_Reduce(&procSum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if (!rank) {
		printf("Answer = %d \n", sum);
	}
	MPI_Finalize();
}



void task8() {
	const int N = 20;
	int a[N], rank, size;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
	int portion = N / size;
	int* ar;
	ar = new int[portion];
	if (!rank)
	{
		for (int i = 0; i < N; i++)
		{
			a[i] = rand() % 40;
			printf("%d ", a[i]);
			if (i == (N - 1)) { printf("\n"); }
		}
		
		printf("Rank : %d ", rank);
		for (int i = 0; i < portion; i++) {
			ar[i] = a[i];
			printf("%d ", ar[i]);
		}
		printf("\n");
		for (int i = 1; i < size; i++) {
			MPI_Send(a+portion*i, portion, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

	}
	else {
		MPI_Recv(ar, portion, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		printf("Rank : %d ", rank);
		for (int i = 0; i < portion; i++) {
			printf("%d ", ar[i]);
		}
		printf("\n");
	}

	printf("\n");
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 1) 
	{
		int *arR, *ans;
		arR = new int[portion];
		ans = new int[size];
		for (int i = 0; i < size; i++) {
			if (i != 1) {
				MPI_Recv(arR, portion, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			}
			else { arR = ar; }
			int startI = i * portion;
			for (int j = 0; j < portion; j++) {
				ans[startI] = arR[j];
				startI++;
			}
		}
		printf("Final: ");
		for (int i = 0; i < N; i++) {
			printf("%d ", ans[i]);
		}
	}
	else {
		MPI_Send(ar, portion, MPI_INT, 1, 1, MPI_COMM_WORLD);
	}


	MPI_Finalize();
}

void task9() {

	
	int rank, size;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	const int N = 2*(size-1);


	MPI_Status status;
	int portion = N / (size - 1);
	int* arR;
	arR = new int[portion];
	int* ans;
	ans = new int[N];
	if (!rank)
	{
		int* a;
		a = new int[N];
		for (int i = 0; i < N; i++)
		{
			a[i] = rand() % 40;
			printf("%d ", a[i]);
		}
		printf("\n");
		for (int i = 1; i < size; i++) {
			MPI_Send(a + portion * (i-1), portion, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		int ind = 0;
		for (int i = size-1; i >0; i--) {
			MPI_Recv(arR, portion, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			for (int j = 0; j < portion; j++) {
				ans[ind] = arR[j];
				ind++;
			}
		}
		printf("Final: ");
		for (int i = 0; i < N; i++) {
			printf("%d ", ans[i]);
		}
	}
	else {
		int portion =N/(size-1);
		int* ar;
		ar = new int[portion];
		int* rev;
		rev = new int[portion];
		MPI_Recv(ar, portion, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		
		for (int i = 0; i < portion; i++) {
			rev[i] = ar[portion - i - 1];
		}
		MPI_Send(rev, portion, MPI_INT, 0, 1, MPI_COMM_WORLD);
		
	}

	MPI_Finalize();
}



void task10() {
	int rank, size, i;
	MPI_Status status;
	double start, end;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	const int N = 1000;
	int* array;
	int* ans;
	array = new int[N];
	ans = new int[N];
	if (!rank) {
		for (i = 0; i < N; i++) {
			array[i] = rand() % 40;
		}
		MPI_Send(array, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(ans, N, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
	}
	else {
		int* a;
		a = new int[N];
		MPI_Recv(a, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Send(a, N, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}


	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();
	if (!rank) { printf("Runtime SEND = %f\n", end - start); }


	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	if (!rank) {
		for (i = 0; i < N; i++) {
			array[i] = rand() % 40;
		}
		MPI_Ssend(array, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(ans, N, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
	}
	else {
		int* a;
		a = new int[N];
		MPI_Recv(a, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Ssend(a, N, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}


	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();
	if (!rank) {
		printf("Runtime SSEND = %f\n", end - start);
	}

	//BSEND
	MPI_Barrier(MPI_COMM_WORLD);

	start = MPI_Wtime();
	int *a = new int[N];
	
	
	if (!rank) {
		for (i = 0; i < N; i++) {
			array[i] = rand() % 40;
		}
		int message_buffer_size = N * sizeof(int) + MPI_BSEND_OVERHEAD;
		int* message_buffer = (int*)malloc(message_buffer_size);
		MPI_Buffer_attach(message_buffer, message_buffer_size);
		MPI_Bsend(array, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(ans, N, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
	}
	else {
		MPI_Recv(a, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		int message_buffer_size = N * sizeof(int) + MPI_BSEND_OVERHEAD;
		int* message_buffer = (int*)malloc(message_buffer_size);
		MPI_Buffer_attach(message_buffer, message_buffer_size);
		MPI_Bsend(a, N, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();
	if (!rank) {
		printf("Runtime BSEND = %f\n", end - start);
	}
	


	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	if (!rank) {
		for (i = 0; i < N; i++) {
			array[i] = rand() % 40;
		}
		MPI_Rsend(array, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(ans, N, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
	}
	else {
		int* a;
		a = new int[N];
		MPI_Recv(a, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Rsend(a, N, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}


	MPI_Barrier(MPI_COMM_WORLD);
	end = MPI_Wtime();
	if (!rank) {
		printf("Runtime RSEND = %f\n", end - start);
	}

	MPI_Finalize();
}

void task11() {;
	int *ans, rank, size;
	int a[4];
	ans = new int[4];
	for (int i = 0; i < 4; i++) {
		a[i] = 0;
	}

	MPI_Status status;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (!rank) {
		MPI_Send(a, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	
		MPI_Recv(ans, 1, MPI_INT, (size-1), 0, MPI_COMM_WORLD, &status);
		ans[0]++;
		printf("%d", ans[0]);
	}
	else {
		
		MPI_Recv(ans, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
		ans[0]++;
		if (rank== (size-1)){ MPI_Send(ans, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); }
		else { MPI_Send(ans, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD); }
	}
	MPI_Finalize();
}

int main()
{
    return EXIT_SUCCESS;
}
