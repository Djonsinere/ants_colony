#include <cstdio>       
#include <vector>       
#include <random>       
#include <algorithm>    
#include <cstddef>      
#include <cmath>        


// Матрица, которая используется для хранения феромонов 
template<typename T>
class Matrix {
public:
    // Конструктор по умолчанию — пустая матрица
    Matrix() : rows_(0), cols_(0), data_() {} 

    // Квадратная матрица размером size x size, заполненная по умолчанию 0
    Matrix(size_t size)
        : rows_(size), cols_(size), data_(size, std::vector<T>(size, T())) {}

    // Для записи
    T& operator()(size_t row, size_t col) {
        return data_[row][col];
    }

    // Для чтения
    const T& operator()(size_t row, size_t col) const {
        return data_[row][col];
    }

    size_t getRows() const { return rows_; }
    size_t getCols() const { return cols_; }

private:
    size_t rows_, cols_;                        // Размеры матрицы
    std::vector<std::vector<T>> data_;          // Данные матрицы
};

// Граф, который используется для хранения расстояний между городами 
template<typename T>
class Graph {
public:
    // Конструктор, создающий граф с size вершинами (матрица смежности)
    Graph(size_t size)
        : size_(size), data_(size, std::vector<T>(size, T{})) {}

    // Для записи
    T& operator()(size_t row, size_t col) {
        return data_[row][col];
    }

    // Для чтение
    const T& operator()(size_t row, size_t col) const {
        return data_[row][col];
    }

    size_t getVertexesCount() const { return size_; }
    
    // Возвращает вес графа, обычно он всегда 1, но можно изменять
    double getGraphWeight() const { return 1.0; }

    // Проверка на пустоту
    bool IsEmpty() const { return size_ == 0; }

private:
    size_t size_;                                // Количество вершин
    std::vector<std::vector<T>> data_;           // Матрица смежности
};

//  Путь муравья, маршрут и длина 
struct AntPath {
    std::vector<size_t> vertices; // Последователность городов
    double distance = 0;          // Длинна пути в целом
};

// Структура муравья 
struct Ant {
    Ant(size_t start)
        : start_location(start), current_location(start), can_continue(true) {}

    size_t start_location;              // Начало
    size_t current_location;            // Текущая точка 
    bool can_continue;                  // Может ли двигаться дальше
    AntPath path;                       // Путь, по которому идёт муравей
    std::vector<size_t> visited;        // Уже посещённые вершины

    // Генерация случайного числа от 0 до 1
    double getRandomChoice();

    // Получение соседей текущей вершины
    std::vector<size_t> getNeighborVertexes(const Graph<double> &g);

    // Совершение одного шага
    void MakeChoice(const Graph<double> &g, const Matrix<double> &p, double alpha, double beta);
};

// Алгоритм муравьиной колонии 
class AntColonyOptimization {
public:
    explicit AntColonyOptimization(const Graph<double> &graph);

    // Решение задачи коммивояжёра
    AntPath SolveSalesmansProblem();

private:
    const double kAlpha_ = 1.0;          // Влияние феромона
    const double kBeta_ = 2.0;           // Влияние расстояния
    const double kPheromone0_ = 1.0;     // Начальное значение феромона
    double kQ_ = 5.0;                    // Константа для обновления феромона
    const double kEvaporation_ = 0.2;    // Коэффициент испарения феромона

    Graph<double> graph_;                // Исходный граф
    Matrix<double> pheromone_;          // Матрица феромонов

    std::vector<Ant> ants_;             // Все муравьи

    void CreateAnts();                  // Создание муравьёв
    void UpdateGlobalPheromone(const Matrix<double> &local_pheromone_update); // Обновляем ферромоны
};

//  Конструктор алгоритма муравьиной колонии 
AntColonyOptimization::AntColonyOptimization(const Graph<double> &graph)
        : graph_(graph), pheromone_(graph.getVertexesCount()), kQ_(0.015 * graph.getGraphWeight()) {

    const size_t kVertexesCount = graph_.getVertexesCount();

    // Инициализация феромонов
    for (size_t row = 0; row != kVertexesCount; ++row)
        for (size_t col = 0; col != kVertexesCount; ++col)
            if (row != col) pheromone_(row, col) = kPheromone0_;
}

//  Создание муравьёв, по одному на каждую вершину 
void AntColonyOptimization::CreateAnts() {
    const auto kAntsCount = graph_.getVertexesCount();
    ants_.clear();
    for (size_t i = 0; i < kAntsCount; ++i)
        ants_.emplace_back(i);
}

//  Обновление феромонов после прохода всех муравьёв 
void AntColonyOptimization::UpdateGlobalPheromone(const Matrix<double> &lpu) {
    for (size_t from = 0, size = lpu.getRows(); from != size; ++from) {
        for (size_t to = 0; to != size; ++to) {
            // Испарение и добавление нового феромона
            pheromone_(from, to) = (1 - kEvaporation_) * pheromone_(from, to) + lpu(from, to);
            if (pheromone_(from, to) < 0.01 && from != to)
                pheromone_(from, to) = 0.01; // Устанавливается минимальное значение ферромона
        }
    }
}

// Получение случайного числа от 0 до 1 
double Ant::getRandomChoice() {
    std::uniform_real_distribution<> dist(0.0, 1.0);
    static std::default_random_engine re(std::random_device{}());
    return dist(re);
}

//  Получение доступных соседей 
std::vector<size_t> Ant::getNeighborVertexes(const Graph<double> &g) {
    std::vector<size_t> vertexes;
    for (size_t to = 0; to != g.getVertexesCount(); ++to) {
        bool edge_is_exist = g(current_location, to) != 0;
        bool vertex_is_unvisited = std::find(visited.begin(), visited.end(), to) == visited.end();
        if (edge_is_exist && vertex_is_unvisited)
            vertexes.push_back(to);
    }
    return vertexes;
}

//  Выбор следующей вершины 
void Ant::MakeChoice(const Graph<double> &g, const Matrix<double> &p, double alpha, double beta) {
    if (path.vertices.empty()) {
        path.vertices.push_back(current_location);
        visited.push_back(current_location);
    }

    std::vector<size_t> neighbor_vertexes = getNeighborVertexes(g);

    // Если нет соседей тогда возвращаемся в начало
    if (neighbor_vertexes.empty()) {
        can_continue = false;
        if (g(current_location, start_location) != 0) {
            path.vertices.push_back(start_location);
            path.distance += g(current_location, start_location);
        }
        return;
    }

    std::vector<double> wish;
    double summary_wish = 0.0f;

    // Расчёт желания муравья пойти в каждый соседний город
    for (auto v : neighbor_vertexes) {
        double t = p(current_location, v); // феромон
        double w = g(current_location, v); // расстояние
        double n = 1 / w;                  // привлекательность
        wish.push_back(std::pow(t, alpha) * std::pow(n, beta));
        summary_wish += wish.back();
    }

    // Расчёт вероятностей перехода
    std::vector<double> choosing_probability(neighbor_vertexes.size());
    std::vector<double> probability;
    for (size_t neighbor = 0; neighbor != neighbor_vertexes.size(); ++neighbor) {
        probability.push_back(wish[neighbor] / summary_wish);
        if (neighbor == 0)
            choosing_probability[neighbor] = probability.back();
        else
            choosing_probability[neighbor] = choosing_probability[neighbor - 1] + probability.back();
    }

    // Выбор следующей вершины
    size_t next_vertex = 0;
    double choose = getRandomChoice();
    for (size_t n = 0; n != neighbor_vertexes.size(); ++n) {
        if (choose <= choosing_probability[n]) {
            next_vertex = neighbor_vertexes[n];
            break;
        }
    }

    path.vertices.push_back(next_vertex);
    path.distance += g(current_location, next_vertex);
    visited.push_back(next_vertex);
    current_location = next_vertex;
}

// Реализация муравьиного алгоритма для TSP 
AntPath AntColonyOptimization::SolveSalesmansProblem() {
    if (graph_.IsEmpty())
        return {};

    constexpr size_t kMaxIterationsWithoutImprovement = 100;
    const size_t kVertexesCount = graph_.getVertexesCount();
    size_t counter = 0;

    AntPath path;
    path.distance = std::numeric_limits<double>::infinity(); // Изначально бесконечная длина

    while (counter++ != kMaxIterationsWithoutImprovement) {
        Matrix<double> local_pheromone_update(kVertexesCount);
        for (size_t i = 0; i < kVertexesCount; ++i)
            for (size_t j = 0; j < kVertexesCount; ++j)
                local_pheromone_update(i, j) = 0.0;

        CreateAnts(); // создаём муравьёв

        for (auto &ant : ants_) {
            while (ant.can_continue)
                ant.MakeChoice(graph_, pheromone_, kAlpha_, kBeta_);

            auto ant_path = ant.path;
            if (ant_path.vertices.size() == kVertexesCount + 1) {
                if (path.distance > ant.path.distance) {
                    path = ant.path;
                    counter = 0;
                }

                for (size_t v = 0; v != ant_path.vertices.size() - 1; ++v)
                    local_pheromone_update(ant_path.vertices[v], ant_path.vertices[v + 1]) += kQ_ / ant_path.distance;
            }
        }

        UpdateGlobalPheromone(local_pheromone_update);
    }
    return path;
}

//  Точка входа в программу 
int main() {

    Graph<double> g(5); // Задание графа из 5 вершин

    // Задаём матрицу расстояний (симметричный граф)
    double dist[5][5] = {
        {0, 2, 9, 10, 7},
        {2, 0, 6, 4, 3},
        {9, 6, 0, 8, 5},
        {10, 4, 8, 0, 6},
        {7, 3, 5, 6, 0}
    };

    // Заполнение графа
    for (size_t i = 0; i < 5; ++i)
        for (size_t j = 0; j < 5; ++j)
            g(i, j) = dist[i][j];

    AntColonyOptimization aco(g);                // Создаём оптимизатор
    AntPath best = aco.SolveSalesmansProblem();  // Запускаем алгоритм

    // Выводим лучший найденный путь
    printf("Best path distance: %f\n", best.distance);
    printf("Path: ");
    for (size_t v : best.vertices) {
        printf("%zu ", v);
    }
    printf("\n");

    return 0;
}