
use gnuplot::*;

pub fn figure_1(low_res_time: &[f64], low_res_brownian: &[f64],
                high_res_time: &[f64], high_res_brownian: &[f64]) {

    let col1 = "#000000";
    let col2 = "#2CA02C";

    let mut fg = Figure::new();
    let ax = fg.axes2d();

    ax.set_x_grid(true);
    ax.set_y_grid(true);

    ax.lines(low_res_time, low_res_brownian, &[Color(col1)]);
    ax.points(low_res_time, low_res_brownian, &[Color(col1), PointSymbol('o')]);


    ax.lines(high_res_time, high_res_brownian, &[Color(col2)]);
    ax.points(high_res_time, high_res_brownian, &[Color(col2), PointSymbol('o')]);



    // set smoother gnuplot terminal and show it
    fg.set_terminal("wxt", "");

    // spawn plot on new thread to let Rust code continue
    std::thread::spawn(move || {
        fg.show().unwrap();
    });
}



enum Color {
    Black,
    Red,
    Green,
    Blue,
    Yellow,
    Cyan,
    Magenta,
    Silver,
}


impl Color {
    fn to_hex(&self) -> String {
        match self {
            Color::Black => String::from("#000000"),
            Color::Red => String::from("#FF0000"),
            Color::Green => String::from("#2CA02C"),
            Color::Blue => String::from("#0000FF"),
            Color::Yellow => String::from("#FFD700"),
            Color::Cyan => String::from("#00FFFF"),
            Color::Magenta => String::from("#D100D1"),
            Color::Silver => String::from("#C0C0C0"),
        }
    }
}


pub fn figure_2(results: &Vec<(usize, Vec<f64>, Vec<f64>)>) {

    let num_plots = results.len();
    let colors = vec![Color::Black, Color::Red, Color::Green,
                      Color::Blue, Color::Yellow, Color::Cyan,
                      Color::Magenta, Color::Silver];

    let mut fg = Figure::new();
    let ax = fg.axes2d();
    ax.set_x_grid(true);
    ax.set_y_grid(true);

    for i in 0..num_plots {
        let col = &colors[i];

        let (_, time, val) = &results[i];


        ax.lines(time, val, &[Color(&col.to_hex())]);
        ax.points(time, val, &[Color(&col.to_hex()), PointSymbol('o'), PointSize(1.0 + (1.0 / ((i+1) as f64)))]);

    }

    // set smoother gnuplot terminal and show it
    fg.set_terminal("wxt", "");

    // spawn plot on new thread to let Rust code continue
    std::thread::spawn(move || {
        fg.show().unwrap();
    });

}
