#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#define TRUE 1
#define FALSE 0

struct problema {
    int tamanho;
    int** elementos;
};

struct nodo {
    int indice;
    int valor;
};

int debug;

void ler_arquivo(struct problema*, char[20]);
void imprimir_solucao(int, int*);
int calcular_custo(struct problema, int*);
void copiar_solucao(int, int*, int*);
int* construir_solucao(struct problema, float, float);
int* mlp_2opt(struct problema, int*);
int* insercao(struct problema, int*);
int* path_relinking(struct problema, int*, int*);
int localizar_elemento(int*, int, int, int);

void selection_sort(struct nodo *, int);
void linha();


void gerar_vizinho_aleatorio(struct problema, int*, int*);
int* sa(struct problema, double, double, int*, int);
int* grasp(struct problema, int, float, float, int);

void imprimir_matriz(int**, int);

int rnd(int, int);

int* inicializar_solucao(int, int*);


//-----------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    int *solucao;
    int iteracoes;
    float alpha_inicial, alpha_final;
    int utilizar_sa;
    time_t inicio;
    struct problema p;
    
    //inicializando o algoritimo de geração de numeros aleatorios
    //com uma semente estatica
    srand(10);
    
    if(argc == 7) {
        //lendo o arquivo da instância
        ler_arquivo(&p, argv[1]);
        
        iteracoes = atoi(argv[2]);
        
        alpha_inicial = atof(argv[3]);
        alpha_final = atof(argv[4]);
        
        utilizar_sa = atoi(argv[5]);
        debug = atoi(argv[6]);
    } else {
        //-- configurações de teste
        ler_arquivo(&p, "/Users/gleissonassis/Dropbox/Mestrado/Implementações/minima-latencia/grasp-minimalatencia/instancias/30_1_100_1000.txt");
     
        iteracoes = 10000;
        alpha_inicial = 0.1;
        alpha_final = 0.9;
        utilizar_sa = 1;
        debug = 1;
    }
    
    //inicializando a contagem de tempo
    inicio = time(NULL);
    
    
    
    if(debug) {
        printf("Instancia : %s\n\n", argv[1]);
    }
    
    solucao = grasp(p, iteracoes, alpha_inicial, alpha_final, utilizar_sa);
    
    if(debug) {
        printf("Caminho = ");
        imprimir_solucao(p.tamanho + 1, solucao);
        printf("Custo = %d\n", calcular_custo(p, solucao));
        printf("Tempo gasto : %d",  (int)(time(NULL) - inicio));
        
        linha();
    } else {
        printf("%s\n", argv[1]);
        imprimir_solucao(p.tamanho + 1, solucao);
        printf("%d\n%d\n", calcular_custo(p, solucao), (int)(time(NULL) - inicio));
    }
    
    //liberando memória para a matriz representando o problema solução
    free(solucao);
    free(p.elementos);
    
    return 0;
}

int rnd(int min, int max) {
    min = ceil(min);
    max = floor(max);
    
    return floor(rand() % (max + 1 - min)) + min;
}

int* grasp(struct problema p, int iteracoes, float alpha_inicial, float alpha_final, int utilizar_sa) {
    int *solucao_atual = inicializar_solucao(p.tamanho, NULL);
    int *solucao_pr = inicializar_solucao(p.tamanho, NULL);
    int* solucao = inicializar_solucao(p.tamanho, NULL);
    int custo_atual;
    int custo = INT_MAX;
    int iteracao;
    time_t inicio = time(NULL);
    
    do {
        solucao_atual = construir_solucao(p, alpha_inicial, alpha_final);
        
        custo_atual = calcular_custo(p, solucao_atual);
        
        if(utilizar_sa) {
            solucao_atual = sa(p, 0.1, custo_atual, solucao_atual, ceil(p.tamanho * 0.2));
        
            custo_atual = calcular_custo(p, solucao_atual);
        
            if(custo_atual < custo) {
                copiar_solucao(p.tamanho, solucao_atual, solucao);
                custo = custo_atual;
            
                if(debug) {
                    printf("Custo melhorado (Simulated Annealing): %d\n", custo);
                    imprimir_solucao(p.tamanho + 1, solucao);
                    printf("Tempo gasto: %d\n", (int)(time(NULL) - inicio));
                    printf("Iteracao: %d\n\n", iteracao);
                }
            }
        }
        
        for(int i = 0; i < p.tamanho; i++) {
            solucao_atual = mlp_2opt(p, solucao_atual);
            
            custo_atual = calcular_custo(p, solucao_atual);
            
            if(custo_atual < custo) {
                copiar_solucao(p.tamanho, solucao_atual, solucao);
                custo = custo_atual;
                
                if(debug) {
                    printf("Custo melhorado (2-opt (%d)): %d\n", i, custo);
                    imprimir_solucao(p.tamanho + 1, solucao);
                    printf("Tempo gasto: %d\n", (int)(time(NULL) - inicio));
                    printf("Iteracao: %d\n\n", iteracao);
                }
            }
        }
        
        solucao_atual = insercao(p, solucao_atual);
        
        custo_atual = calcular_custo(p, solucao_atual);
        
        if(custo_atual < custo) {
            copiar_solucao(p.tamanho, solucao_atual, solucao);
            custo = custo_atual;
            
            if(debug) {
                printf("Custo melhorado (insercao): %d\n", custo);
                imprimir_solucao(p.tamanho + 1, solucao);
                printf("Tempo gasto: %d\n", (int)(time(NULL) - inicio));
                printf("Iteracao: %d\n\n", iteracao);
                
            }
        }
        
        custo_atual = calcular_custo(p, solucao_atual);
        
        if(custo_atual > custo){
            solucao_pr = path_relinking(p, solucao, solucao_atual);
            custo_atual = calcular_custo(p, solucao_pr);
            
            if(custo_atual < custo) {
                copiar_solucao(p.tamanho, solucao_pr, solucao);
                custo = custo_atual;
                
                if(debug) {
                    printf("Custo melhorado (PATH): %d\n", custo);
                    imprimir_solucao(p.tamanho+1, solucao);
                    printf("Tempo gasto: %d\n", (int)(time(NULL) - inicio));
                    printf("Iteracao: %d\n\n", iteracao);
                    
                }
            }
        }
        
        iteracoes--;
        iteracao++;
    }while(iteracoes > 0);
    
    return solucao;
}

//-----------------------------------------------------------------------------

int* construir_solucao(struct problema p, float percentual_inicial, float percentual_final) {
    int iv, indice_selecionado, indice_selecionado2;
    int *inserido;
    struct nodo *vizinhos;
    int numero_candidatos;
    float taxa_crescimento, percentual_atual;
    int* solucao = inicializar_solucao(p.tamanho, NULL);
    
    //estabelecendo o numero de candidatos a entrar no solucao
    numero_candidatos = ceil(percentual_inicial * p.tamanho);
    
    taxa_crescimento = (percentual_final - percentual_inicial) / p.tamanho;
    
    percentual_atual = percentual_inicial + taxa_crescimento;
    
    //alocando memoria para o array que ira informar se um elemento
    //ja foi inserido no solucao ou nao e inicializando os valores
    inserido = malloc(p.tamanho * sizeof(int));
    for(int i = 0; i < p.tamanho; i++) {
        inserido[i] = FALSE;
    }
    
    //alocando memoria para o array que irá armazenar informações sobre a vizinhança
    //do elemento atual
    vizinhos = (struct nodo*) malloc((p.tamanho) * sizeof(struct nodo));
    
    solucao[0] = 0;
    inserido[0] = TRUE;
    
    for(int i = 0; i < p.tamanho; i++) {
        //indice do vizinho atual;
        iv = 0;
        
        //construindo a lista de vizinhos e seus respectivos valores
        for(int j = 0; j < p.tamanho; j++) {
            //nao pode ser selecionado um elemento que ja esteja no solucao
            if(!inserido[j]) {
                vizinhos[iv].indice = j;
                vizinhos[iv].valor = p.elementos[i][j];
                
                iv++;
            }
        }
        
        if(iv == 0) {
            solucao[i + 1] = 0;
            //printf("Vertice inserido no solucao: %d\n", solucao[i + 1]);
        } else {
            //ordenando a lista de vizinhos por custo da aresta
            selection_sort(vizinhos, iv);
            
            //selecionando um elemento aleatorio para entrar no solucao
            do {
                if(numero_candidatos > iv) {
                    indice_selecionado = rand() % iv;
                    indice_selecionado2 = rand() % iv;
                } else {
                    indice_selecionado = rand() % numero_candidatos;
                    indice_selecionado2 = rand() % numero_candidatos;
                }
                
                if(!inserido[vizinhos[indice_selecionado2].indice] && vizinhos[indice_selecionado].valor > vizinhos[indice_selecionado2].valor) {
                    indice_selecionado = indice_selecionado2;
                }
            } while(inserido[vizinhos[indice_selecionado].indice] == TRUE);
            solucao[i + 1] = vizinhos[indice_selecionado].indice;
            inserido[vizinhos[indice_selecionado].indice] = TRUE;
            
            //printf("Vertice inserido no solucao: %d\n", solucao[i + 1]);
            
            numero_candidatos = ceil(percentual_atual * p.tamanho);
            percentual_atual += taxa_crescimento;
            
        }
    }
    
    //liberando memoria alocada para os arrays
    free(inserido);
    free(vizinhos);
    
    return solucao;
}

/*
 * Function: gerar_vizinho_aleatorio
 * -----------------------------------------------------------------------------
 *   Função que perturba uma solução com o objetivo de exporar de forma melhor
 *   o espaço de soluções.
 *
 *   s: estrutura que representa uma construção de solução para o problema.
 */
void gerar_vizinho_aleatorio(struct problema p, int* solucao, int* solucao_destino) {
    int tmp, i, j;
    
    copiar_solucao(p.tamanho, solucao, solucao_destino);
    
    if(rnd(0, 50) > 25) {
        i = rnd(1, p.tamanho - 1);
        j = rnd(1, p.tamanho - 1);
        
        tmp = solucao_destino[i];
        solucao_destino[i] = solucao_destino[j];
        solucao_destino[j] = tmp;
    } else {
        i = rnd(1, p.tamanho - 2);
        
        tmp = solucao_destino[i];
        solucao_destino[i] = solucao_destino[i + 1];
        solucao_destino[i + 1] = tmp;
    }
}

/*
 * Function: sa
 * -----------------------------------------------------------------------------
 *   Implementação da metaheurística Simulated Annealing.
 *
 *   a: fator de redução de temperatura.
 *   T0: temperatura inicial.
 *   s: estrutura que representa uma construção de solução para o problema.
 *   SAmax: quantidade de perturbações em uma iteração do Simulated Annealing
 */
int* sa(struct problema p, double a, double T0, int* solucao_inicial, int perturbarcoes) {
    int* s;
    int* sl;
    int* s_estrela;
    
    int IterT = 0;
    double T = T0;
    double delta = 0;
    
    
    //alocando memória para o solucao que será criado
    s = inicializar_solucao(p.tamanho, solucao_inicial);
    sl = inicializar_solucao(p.tamanho, solucao_inicial);
    s_estrela = inicializar_solucao(p.tamanho, solucao_inicial);
    
    //copiando para as soluções temporárias a solução inicial passada como
    //parâmetro
    
    
    while(T > 0.0001) {
        while(IterT < perturbarcoes) {
            IterT++;
            
            //A geração de um vizinho aleatório simula a perturbação de uma
            //solução
            gerar_vizinho_aleatorio(p, s, sl);
            
            int fsl = calcular_custo(p, sl);
            int fs = calcular_custo(p, s);
            delta = fsl - fs;
            
            //se a solução tem um custo melhor que a solução corrente ela é
            //automaticamente aceita. Observe que a solução corrente não é
            //necessáriamente a melhor solução (s_estrela)
            if(delta < 0) {
                copiar_solucao(p.tamanho, sl, s);
                
                int fs_estrela = calcular_custo(p, s_estrela);
                
                //Se a solução corrente é melhor de todos até o momento ela
                //passa ser s_estrela
                if(fsl < fs_estrela) {
                    copiar_solucao(p.tamanho, sl, s_estrela);
                }
            } else {
                double x = rand() / RAND_MAX;
                
                //aceitando uma solução de piora dada a função de probabilidade
                if(x < pow(M_E, -delta / T)) {
                    copiar_solucao(p.tamanho, sl, s);
                }
            }
        }
        
        //Reduzindo a temperatura através do fator de redução e zerando a contagem
        //de perturbações
        IterT = 0;
        T = T - (a * T);
    }
    
    copiar_solucao(p.tamanho, s_estrela, s);
    
    free(s);
    free(sl);
    
    return s_estrela;
}

//-----------------------------------------------------------------------------

int localizar_elemento(int* solucao, int tamanho, int valor, int posicao_inicial) {
    int i = posicao_inicial;
    
    for(; i < tamanho; i++) {
        if(solucao[i] == valor) {
            return i;
        }
    }
    
    return -1;
}

//-----------------------------------------------------------------------------

int* path_relinking(struct problema p, int* c1, int* c2) {
    int i, j;
    int custo, tmp_custo;
    int pos;
    int tmp;
    int *solucao_resultado;
    int *solucao_tmp;
    
    
    //alocando memória para as soluções que serão criadas
    solucao_resultado = inicializar_solucao(p.tamanho, c1);
    solucao_tmp = inicializar_solucao(p.tamanho, c2);
    
    solucao_tmp = malloc((p.tamanho + 1) * sizeof(int));
    
    copiar_solucao(p.tamanho, c1, solucao_resultado);
    copiar_solucao(p.tamanho, c2, solucao_tmp);
    
    custo = calcular_custo(p, c1);
    
    for(i = 1; i < p.tamanho; i++) {
        if(c1[i] == solucao_tmp[i]) {
            continue;
        }
        
        pos = localizar_elemento(solucao_tmp, p.tamanho, c1[i], i + 1);
        
        for(j = pos; j > i; j--) {
            tmp = solucao_tmp[j];
            solucao_tmp[j] = solucao_tmp[j - 1];
            solucao_tmp[j - 1] = tmp;
            
            tmp_custo = calcular_custo(p, solucao_tmp);
            
            if(tmp_custo < custo) {
                custo = tmp_custo;
                copiar_solucao(p.tamanho, solucao_tmp, solucao_resultado);
            }
        }
    }
    
    free(solucao_tmp);
    
    return solucao_resultado;
}

//-----------------------------------------------------------------------------

int* insercao(struct problema p, int* solucao) {
    int custo, custo_insercao;
    int* solucao_resultado;
    int * solucao_tmp;
    int tmp;
    
    custo = calcular_custo(p, solucao);
    
    
    solucao_resultado = malloc((p.tamanho + 1) * sizeof(int));
    solucao_tmp = malloc((p.tamanho + 1) * sizeof(int));
    
    copiar_solucao(p.tamanho, solucao, solucao_resultado);
    copiar_solucao(p.tamanho, solucao, solucao_tmp);
    
    for(int i = 1; i < p.tamanho; i++) {
        copiar_solucao(p.tamanho, solucao, solucao_tmp);
        
        for(int j = i + 1; j < p.tamanho; j++) {
            tmp = solucao_tmp[i];
            solucao_tmp[i] = solucao_tmp[j];
            solucao_tmp[j] = tmp;
            custo_insercao = calcular_custo(p, solucao_tmp);
            
            if(custo_insercao < custo) {
                custo = custo_insercao;
                copiar_solucao(p.tamanho, solucao_tmp, solucao_resultado);
            }
        }
    }
    
    //liberando memória alocada para o solucao temporária
    free(solucao_tmp);
    free(solucao);
    
    return solucao_resultado;
}

//-----------------------------------------------------------------------------

int* mlp_2opt(struct problema p, int* solucao) {
    int i,j;
    int custo, custo_swap;
    int* solucao_swap;
    int * solucao_tmp;
    int tmp;
    
    custo = calcular_custo(p, solucao);
    
    solucao_swap = malloc((p.tamanho + 1) * sizeof(int));
    solucao_tmp = malloc((p.tamanho + 1) * sizeof(int));
    
    copiar_solucao(p.tamanho, solucao, solucao_swap);
    copiar_solucao(p.tamanho, solucao, solucao_tmp);
    
    for(i = 1; i < p.tamanho; i++) {
        copiar_solucao(p.tamanho, solucao, solucao_tmp);
        
        for(j = i + 1; j < p.tamanho; j++) {
            tmp = solucao_tmp[i];
            solucao_tmp[i] = solucao_tmp[j];
            solucao_tmp[j] = tmp;
            custo_swap = calcular_custo(p, solucao_tmp);
            
            if(custo_swap < custo) {
                custo = custo_swap;
                copiar_solucao(p.tamanho, solucao_tmp, solucao_swap);
            }
        }
    }
    
    //liberando memória alocada para o solucao temporário
    free(solucao_tmp);
    free(solucao);
    
    return solucao_swap;
}

//-----------------------------------------------------------------------------

void ler_arquivo(struct problema* p, char arquivo[2000]) {
    FILE* fp;
    int t;
    
    fp = fopen(arquivo, "r");
    fscanf(fp, "%d %d\n\n", &p->tamanho, &t);
    
    //alocando espaço para a matriz de adjacencia
    p->elementos = malloc(p->tamanho * p->tamanho * sizeof(int));
    
    
    //pulando as linhas de 1s
    for(int i = 0; i < p->tamanho; i++) {
        fscanf(fp, "%d ", &t);
    }
    
    //percorrendo os elementos da matriz de ajdacencia que estao no arquivo
    for(int i = 0; i < p->tamanho; i++) {
        p->elementos[i] = malloc(p->tamanho * sizeof(int));
        
        for(int j = 0; j < p->tamanho; j++) {
            p->elementos[i][j] = 0;
            fscanf(fp, "%d ", &p->elementos[i][j]);
        }
    }
}

//-----------------------------------------------------------------------------

void selection_sort(struct nodo *array, int n) {
    int i, j;
    int min;
    struct nodo temp;
    
    for(i = 0; i < n - 1; i++) {
        min=i;
        for(j = i + 1; j < n; j++) {
            if(array[j].valor < array[min].valor) {
                min = j;
            }
        }
        
        temp = array[i];
        array[i] = array[min];
        array[min] = temp;
    }
}

//-----------------------------------------------------------------------------

void copiar_solucao(int n, int* origem, int* destino) {
    int i;
    
    for(i = 0; i < n + 1; i++) {
        destino[i] = origem[i];
    }
}

//-----------------------------------------------------------------------------

int calcular_custo(struct problema p, int* solucao) {
    int custo = 0;
    int i;
    
    for(i = 0; i < p.tamanho; i++) {
        custo += p.elementos[solucao[i]][solucao[i + 1]] * (p.tamanho - i);
    }
    
    return custo;
}

//-----------------------------------------------------------------------------

void imprimir_solucao(int n, int* solucao) {
    int i;
    
    for(i = 0; i < n; i++) {
        printf("%d ", solucao[i]);
    }
    printf("\n");
}

//-----------------------------------------------------------------------------

void linha() {
    int i;
    printf("\n");
    for(i = 0; i < 80; i++) printf("_");
    printf("\n");
}

void imprimir_matriz(int** elementos, int n) {
    linha();
    printf("Matriz\n\n");
    
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%d ", elementos[i][j]);
        }
        
        printf("\n");
    }
    
    linha();
}


int* inicializar_solucao(int n, int* solucao_base) {
    int* nova_solucao;
    
    nova_solucao = malloc((n + 1) * sizeof(int));
    
    if(solucao_base) {
        copiar_solucao(n, solucao_base, nova_solucao);
    }
    
    return nova_solucao;

}
