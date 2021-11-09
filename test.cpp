#include <stdio.h>
#include <omp.h>
#include <mpi.h>

#define MASTER 0
#define START 1
#define QUERY 2
#define ANSWER 3
#define FINISH 4

int main(int argc, char* argv[]) {
    int num_tasks, task_id, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int iam = 0, np = 1;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Get_processor_name(processor_name, &name_len);

    if (task_id == 0) {
        int tmp;
        MPI_Status status;
        printf("Thread %d on master processor %s, rank %d out of %d processors\n", iam, processor_name, task_id, num_tasks);
        for (int worker = 1; worker < num_tasks; worker++) {
            MPI_Send(&task_id, 1, MPI_INT, worker, START, MPI_COMM_WORLD);
        }
        for (int worker = 1; worker < num_tasks; worker++) {
            MPI_Recv(&tmp, 1, MPI_INT, worker, FINISH, MPI_COMM_WORLD, &status);
        }
        for (int worker = 1; worker < num_tasks; worker++) {
            MPI_Send(&task_id, 1, MPI_INT, worker, QUERY, MPI_COMM_WORLD);
        }
        for (int worker = 1; worker < num_tasks; worker++) {
            MPI_Recv(&tmp, 1, MPI_INT, worker, FINISH, MPI_COMM_WORLD, &status);
        }
    }
    else {
        int start;
        MPI_Status status;
        MPI_Recv(&start, 1, MPI_INT, MASTER, START, MPI_COMM_WORLD, &status);

        int a[4][2] = {{task_id,task_id},{task_id,task_id},{task_id,task_id},{task_id,task_id}};
        #pragma omp parallel default(shared) private(iam, np) num_threads(2)
        {
            np = omp_get_num_threads();
            iam = omp_get_thread_num();
            int tmp = 0;

            // Calculation thread
            if (iam == 0) {
                MPI_Status t0_status;
                for (int worker = 1; worker < num_tasks; worker++) {
                    if (worker == task_id) continue;
                    for (int j = 0; j < 2; j++) {
                        MPI_Send(&task_id, 1, MPI_INT, worker, QUERY, MPI_COMM_WORLD);
                        MPI_Recv(&tmp, 1, MPI_INT, worker, ANSWER, MPI_COMM_WORLD, &t0_status);
                        a[worker-1][j] = tmp;
                    }
                }
                printf("Thread %d on worker processor %s, rank %d out of %d processors, a = {{%d,%d},{%d,%d},{%d,%d},{%d,%d}}\n", iam, processor_name, task_id, num_tasks, a[0][0], a[0][1], a[1][0],a[1][1], a[2][0], a[2][1], a[3][0], a[3][1]);
                int finish = 1;
                MPI_Send(&finish, 1, MPI_INT, MASTER, FINISH, MPI_COMM_WORLD);
            }

            // Listenner thread
            else {
                MPI_Status t1_status;
                while (true) {
                    MPI_Recv(&tmp, 1, MPI_INT, MPI_ANY_SOURCE, QUERY, MPI_COMM_WORLD, &t1_status);
                    if (t1_status.MPI_SOURCE == MASTER) break;
                    MPI_Send(&task_id, 1, MPI_INT, t1_status.MPI_SOURCE, ANSWER, MPI_COMM_WORLD);
                }
                int finish = 1;
                MPI_Send(&finish, 1, MPI_INT, MASTER, FINISH, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Finalize();
}