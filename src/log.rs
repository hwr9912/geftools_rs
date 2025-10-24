use chrono::Local;

pub fn log_msg(msg: &str) {
    let now = Local::now();
    eprintln!("[{}] {}", now.format("%Y-%m-%d %H:%M:%S%.3f %:z"), msg);
}
