import React, { useState } from "react";
import Head from "./components/header";
import { BrowserRouter as Router, Route, Routes } from "react-router-dom";
import MoleVisualizer from "./components/moleVisualizer";
import Instructions from "./components/instructions";
import "./App.css";

function App() {
  
  return (
   
      <Router>
        <Routes>
        <Route path="/" element={<MoleVisualizer />} />
        <Route path="/Head" element={<Head />} />
        <Route path="/Instructions" element={<Instructions />} />

        </Routes>
      </Router>

  );
}

export default App;
